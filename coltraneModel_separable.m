%function OUT = coltraneModel(forcing0,p);

% out = coltraneModel(forcing0);
% ... = coltraneModel(forcing0,p);
%
% time-integrates one population with a single set of traits, but variable
% spawning date tspawn.
%
% For original model description, see Banas et al, Front. Mar. Res., 2016.
% (http://neilbanas.com/projects/coltrane)
%
% this version is an apporximate solution which doesn't require iteration across
% all state variables in time--an exercise in making the model cleanly 
% hierarchical, so that one can derive predictions about development alone (D),
% then size and time evolution in surviving cohorts (D,W), and then mortality, 
% survivorship, and population dynamics (D,W,N). This will also make it possible
% for N to be density-dependent.


%if nargin<2, p = coltraneParams; end
%if nargin<1, forcing0 = coltraneForcing('simple-g'); end

v.tspawn = 0 : p.dt_tspawn : 365; % initial spawn date
NC = length(v.tspawn);

% construct forcing for the whole cohort, interpolating and staggering in time
t = 0:p.dt:((p.Nyears+1)*365);
N = length(t);
forcing.t = repmat(t(:),[1 NC]);
tc = mod(forcing.t,365);
forcing.T =  interp1(forcing0.t, forcing0.T,  tc);
forcing.Tb = interp1(forcing0.t, forcing0.Tb, tc);
forcing.P =  interp1(forcing0.t, forcing0.P,  tc);

% initialize storage
v.lnN = nan.*ones(N,NC); % log number of individuals per starting individual
v.D = nan.*ones(N,NC); % developmental state, 0..1
v.W = nan.*ones(N,NC); % structural biomass
v.a = nan.*ones(N,NC); % activity level
v.m = zeros(N,NC); % mortality rate

% pick a guess at adult size based on mean temperature in the forcing,
% and pick a corresponding guess at egg size based on this. These used to be
% called p.Wa0 and p.We0.
T_nominal = mean(forcing.T(:,1));
qoverq = (p.Q10g./p.Q10d).^(T_nominal./10);
co = (1-p.theta) .* p.GGE_nominal .* (1-p.Df) .* qoverq;
v.Wa_theo = (co .* p.I0 ./ p.u0) .^ (1./(1-p.theta));
v.We_theo = p.r_ea .* v.Wa_theo.^p.exp_ea;

% initial conditions
for i=1:NC
	f = find(t == v.tspawn(i));	
	v.lnN(f,i) = 0;
	v.D(f,i) = 0;
	v.W(f,i) = v.We_theo;
	v.a(f,i) = 1;
end


% activity: a(t)  ------------------------------------
sat = forcing.P ./ (p.Ks + forcing.P);
C = 1;
sat_crit = p.rm .* (1-p.rb) ./ p.r_assim + C .* p.m0./p.I0./p.r_assim;
v.a = double(sat >= sat_crit);


% temperature response factors -----------------------
Teff = v.a .* forcing.T + (1-v.a) .* forcing.Tb;
qd = (p.Q10d^0.1) .^ Teff;
qg = (p.Q10g^0.1) .^ Teff;


% development: D(t) ----------------------------------
isalive = forcing.t >= repmat(v.tspawn,[N 1]); % spawned yet?
dDdt = isalive .* p.u0 .* qd; % nonfeeding formula
v.D = cumsum(dDdt) .* p.dt;
isfeeding = v.D >= p.Df;
dDdt_feeding = isalive .* p.u0 .* qd .* v.a .* sat;
dDdt(isfeeding) = dDdt_feeding(isfeeding); % full formula 
v.D = cumsum(dDdt) .* p.dt;
v.D(v.D>1) = 1;
isfeeding = v.D >= params.Df;
	% because cumsum() is one timestep off from a true forward integration
% development time
isadult = v.D >= 1;
ta = v.t;
ta(~isadult) = inf;
tj = v.t;
tj(~(v.D>0)) = inf;
v.t0eff = min(tj); % effective spawning date, discounting unrealistic long waits
v.dtdev = min(ta) - v.t0eff;
	% development time, from when D starts advancing beyond 0 (t0eff) to D=1
v.dtdev(~any(isadult)) = nan;


% juvenile growth: W(t) ------------------------------
% approximate curve of W(D), valid for a=1
% this is used for the allometry in the calculation of net gain
isgrowing = isfeeding & v.D < 1;
satmean = sum(sat .* isgrowing) ./ sum(isgrowing);
c = (1-p.theta) .* qg ./ qd .* p.I0 ./ p.u0 .* (p.r_assim - p.rm ./ satmean);
Wapprox = v.W;
Wapprox(isgrowing) = ((v.We_theo^(1-p.theta)) + ...
	c(isgrowing) .* (v.D(isgrowing) - p.Df)) .^ (1/(1-p.theta));
% net gain over development (here stored as gain * W)
ImaxW = qg .* p.I0 .* Wapprox .^ p.theta;
	% only used at (isgrowing=1), not necessarily accurate beyond that
astar = p.rb + (1-p.rb) .* v.a;
GW = zeros(size(v.W));
GW(isgrowing) = ImaxW(isgrowing) .* ...
	(v.a(isgrowing) .* p.r_assim .* sat(isgrowing) - p.rm .* astar(isgrowing));
% integrate to get the correct W over development
v.W = v.We_theo + cumsum(GW) .* p.dt;
v.W = max(0,v.W);
% adult size Wa
last = v.D(1:end-1,:) < 1 & v.D(2:end,:)==1;
v.Wa = max(v.W(1:end-1,:) .* last);
isadult = v.D >= 1;
Wa2 = repmat(v.Wa,[N 1]);
v.W(isadult) = Wa2(isadult);
% calculate net gain during adulthood (used later to calculate egg production)
% and fill in juvenile ingestion, metabolism, and net gain while we're at it
Imax = qg .* p.I0 .* v.W.^(p.theta-1);
Imax(~isfeeding) = 0;
Imax(v.W<v.We_theo) = 0;
v.I = v.a .* p.r_assim .* sat .* Imax;
	% this was r_assim * I, as opposed to I, in Coltrane 1.0
v.M = p.rm .* astar .* Imax;
v.G = I - M;


% mortality and survivorship: N(t) ------------------
v.m = p.m0 .* qg .* v.a .* v.W.^(p.theta-1); % mort. rate at T, size
v.m(~isalive) = 0;
v.lnN = cumsum(-v.m) .* p.dt;
% adult recruitment: survivorship at D=1
lnNa = v.lnN;
lnNa(~isadult) = nan;
v.Na = exp(max(lnNa));


% egg production and fitness ------------------------
% one-generation calculation; ignores timing and internal life-history mismatch
Einc = C.G .* C.W; % energy put into income egg prod.
Einc(~isadult) = 0;
C.Finc = sum(Einc.*exp(C.lnN)) .* params.dt ./ C.We_theo;
	% lifetime eggs per egg, from adult accumulation of energy (i.e., income)
C.Fcap = 0.5 .* C.Wa .* C.Na ./ C.We_theo;
	% lifetime eggs per egg, from juvenile accumulation of energy, stored
	% in the form of body mass W (i.e., capital)
	% note the magic number 0.5!
C.F = C.Finc + C.Fcap;
C.capfrac = C.Fcap ./ C.F;


% organise output into structures ----------------------------------------------
OUT.p = p;
OUT.forcing = forcing;
OUT.var = v;


% postprocessing ---------------------------------------------------------------
%OUT.pot = coltranePostproc_phi(v,p,forcing);
%[OUT.cohort, OUT.pop, OUT.routine] = coltranePostproc_ER(v,p,forcing);

