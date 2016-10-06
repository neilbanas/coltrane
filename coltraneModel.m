function [v,p,forcing] = coltraneModel(forcing0,p);

% [v,p,forcing] = coltraneModel(forcing0);
%           ... = coltraneModel(forcing0,p);
%
% time-integrates one population with a single set of traits, but variable
% spawning date tspawn. Calculates summary metrics for both the phi and ER
% versions of the model.
%
% For full model description, see Banas et al (2016) Front. Mar. Res., submitted
% (http://neilbanas.com/projects/coltrane)


if nargin<2, p = coltraneParams; end
if nargin<1, forcing0 = coltraneForcing('simple-g'); end

p.tspawn = 0 : p.dt_tspawn : 365; % initial spawn date
NC = length(p.tspawn);

% construct forcing for the whole cohort, interpolating and staggering in time
t = 0:p.dt:((p.Nyears+1)*365);
N = length(t);
forcing.t = repmat(t(:),[1 NC]);
forcing.tc = mod(forcing.t,365);
forcing.T =  interp1(forcing0.t, forcing0.T,  forcing.tc);
forcing.Tb = interp1(forcing0.t, forcing0.Tb, forcing.tc);
forcing.P =  interp1(forcing0.t, forcing0.P,  forcing.tc);

% initialize storage
v.lnN = nan.*ones(N,NC); % log number of individuals per starting individual
v.D = nan.*ones(N,NC); % developmental state, 0..1
v.S = nan.*ones(N,NC); % structural biomass
v.phi = nan.*ones(N,NC); % potential reserves/eggs
v.a = nan.*ones(N,NC); % activity level
v.R = nan.*ones(N,NC); % reserve biomass
v.Einc = zeros(N,NC); % income egg production
v.Ecap = zeros(N,NC); % capital egg production
v.m = zeros(N,NC); % mortality rate
v.Fin = zeros(N,NC); % ingestion flux
v.Fmet = zeros(N,NC); % metabolism flux

% initial conditions
for i=1:NC
	f = find(t == p.tspawn(i));	
	v.lnN(f,i) = 0;
	v.D(f,i) = 0;
	v.S(f,i) = p.We0;
	v.phi(f,i) = 0;
	v.R(f,i) = 0;
	v.a(f,i) = 1;
end


% ------------------------------------------------------------------------------
% main loop --------------------------------------------------------------------
small = 1e-6;
for n=1:N-1
	sp = forcing.t(n,:) >= p.tspawn; % the ones that have been spawned so far
	% food-saturation and temperature-response coefficients
	sat = forcing.P(n,:) ./ (p.Ks + forcing.P(n,:));
	Teff = v.a(n,:) .* forcing.T(n,:) + (1-v.a(n,:)) .* forcing.Tb(n,:);
	qd = (p.Q10d^0.1) .^ Teff;
	qg = (p.Q10g^0.1) .^ Teff;
	Sthm1 = max(p.We0, v.S(n,:)) .^ (p.theta-1);
	% development
	u = p.u0 .* qd .* sat .* v.a(n,:);
	nf = v.D(n,:) < p.Df; % flag for non-feeding stage
	u(nf) = p.u0 .* qd(nf);
	v.D(n+1,sp) = v.D(n,sp) + u(sp).*p.dt;
	v.D(n+1,v.D(n+1,:)>1) = 1;
	% assimilation
	Imax = p.I0 .* qg .* Sthm1; % max ingestion rate at T, size
	Fin = p.r_assim .* Imax .* sat .* v.a(n,:);
	Fin(nf) = 0; % no gain in the non-feeding stages
	% metabolism and net gain
	Fmet = p.rm .* Imax .* (p.rb + (1-p.rb).*v.a(n,:));
	Fmet(nf) = 0; % G = 0 for non-feeding stages
	G = Fin - Fmet;	
	GSdt = G .* v.S(n,:) .* p.dt;
	% allocation of positive net gain to S and phi
	fs = (min(v.D(n,:),1)-1) ./ (p.Ds(1,:) - 1 - small);
	fs = min(1, max(0, fs.^p.alpha_s)); % fs=1 means entirely to/from S
	GSdtpos = max(0,GSdt);
	v.S(n+1,sp) = v.S(n,sp) + fs(sp) .* GSdtpos(sp);
	v.phi(n+1,sp) = v.phi(n,sp) + (1-fs(sp)) .* GSdtpos(sp);
	% negative net gain is taken from phi
	GSdtneg = min(0,GSdt);
	v.phi(n+1,sp) = v.phi(n+1,sp) + GSdtneg(sp);	
	% determine diapause (instantaneous assessment; fluctuates step by step)
	C = max(0, 1 + min(p.r_phi_max, v.phi(n,:))./v.S(n,:));
	sat_crit = p.rm .* (1-p.rb) ./ p.r_assim + C .* p.m0./p.I0./p.r_assim;
	dia = v.D(n,:) > p.Ddia & sat < sat_crit;
	v.a(n+1,sp & ~dia) = 1;
	v.a(n+1,sp & dia) = 0; % assume 0 activity during diapause
	% mortality
	m = p.m0 .* qg .* Sthm1 .* v.a(n,:); % mortality rate at T, size
	v.m(n+1,:) = m;
	v.lnN(n+1,sp) = v.lnN(n,sp) - m(sp).*p.dt;
	
	% now subdivide phi into R, Einc, and Ecap.
	% start by giving R the same increment as phi
	v.R(n+1,sp) = v.R(n,sp) + (v.phi(n+1,sp) - v.phi(n,sp));
	% has egg production started?
	egg = v.D(n,:) >= 1 ...
	    & forcing.t(n,:) >= p.tegg;
	% calculate egg production and subtract from R
	v.Einc(n,sp & egg) = p.finc .* (1-fs(sp & egg)) .* GSdtpos(sp & egg)./p.dt;	
	Emax = p.r_assim .* Imax .* v.S(n,:);
	v.Emax(n+1,:) = Emax;
	v.R(n+1,sp) = v.R(n+1,sp) - v.Einc(n,sp).*p.dt;
	v.Ecap(n,sp & egg) = p.fcap .* min(max(0,v.R(n+1,sp & egg)./p.dt), ...
						           Emax(sp & egg) - v.Einc(n,sp & egg));
	v.R(n+1,sp) = v.R(n+1,sp) - v.Ecap(n,sp).*p.dt;
	
	% save fluxes for use in energy budget
	v.Fin(n,:) = Fin;
	v.Fmet(n,:) = Fmet;
end
v.lnN(isnan(v.lnN)) = -Inf;
% blank everything after starvation limit is reached in the potential model
for i=1:NC
	f = find(v.D(:,i) > p.Df & v.phi(:,i) < -p.r_starv .* v.S(:,i));
	if ~isempty(f)
		v.D(f(1):end,i) = nan;
		v.S(f(1):end,i) = nan;
		v.phi(f(1):end,i) = nan;
		v.lnN(f(1):end,i) = -Inf;
		v.R(f(1):end,i) = nan;
		v.Einc(f(1):end,i) = nan;
		v.Ecap(f(1):end,i) = nan;
	end
end
% also blank R, Einc, Ecap after starvation limit is reached in the ER model
for i=1:NC
	f = find(v.D(:,i) > p.Df & v.R(:,i) < -p.r_starv .* v.S(:,i));
	if ~isempty(f)
		v.R(f(1):end,i) = nan;
		v.Einc(f(1):end,i) = 0;
		v.Ecap(f(1):end,i) = 0;
	end
end


% ------------------------------------------------------------------------------
% potential model postprocessing -----------------------------------------------

Ny = p.Nyears+1; % Ny = record length, p.Nyears = time allowed for the 
				 % latest-spawned generation to develop
t = forcing.t(:,1);
F0 = v.phi ./ p.We0 .* exp(v.lnN) .* double(v.D >= 1);
	% F0(i,j) is potential eggs ("fitness") at time _forcing.t(i,j)_ for a 
	% cohort spawned on _p.tspawn(j)_
pot.F0 = F0;
pot.t = t;
pot.tspawn = p.tspawn;

% rejigger F0 into a square matrix, covering Nyears+1 x Nyears+1,
% with timebase _tspawn_ on both axes
dt_per_dtspawn = (p.tspawn(2)-p.tspawn(1)) / (t(2)-t(1));
N1y = 365 / (p.tspawn(2)-p.tspawn(1)); % = length(p.tspawn)-1, normally
t1 = t(1:dt_per_dtspawn:end);
F1 = F0(1:dt_per_dtspawn:end,:);
Fn = [];
for i=1:p.Nyears+1
	Fn = [Fn [zeros(N1y*(i-1),N1y); F1(1:end-N1y*(i-1),1:N1y)]];
end
Fn(isnan(Fn)) = 0;
v.pot.Fn = Fn;
v.pot.t1 = t1;

for j=1:length(p.tspawn)
	pot.Fmax(j) = max(F0(:,j)); % max fitness of a 1-generation strategy
	f = find(mod(t-p.tspawn(j),365)==0 & t>p.tspawn(j));
	f = f(1:p.Nyears);
	pot.Fmult(:,j) = F0(f,j); % fitness at integer # years from spawning
	% find optimal sequence of generations
	fs_new = [find(t1==p.tspawn(j)) find(t1==p.tspawn(j)+p.Nyears*365)];
	fs = [];
	while length(fs_new) > length(fs)
		fs = fs_new;
		fs_new = [];
		for k=1:length(fs)-1
			fs_new = [fs_new fs(k)];
			Fkold = Fn(fs(k+1),fs(k));
			% try to insert a new generation in between fs(k) and fs(k+1),
			% see if it increases fitness
			Fknew = [];
			for m = 1 : (fs(k+1)-fs(k)-1)
				Fknew(m) = Fn(fs(k)+m,fs(k)) * Fn(fs(k+1),fs(k)+m);
			end
			if max(Fknew) > Fkold
				mmax = find(Fknew==max(Fknew));
				mmax = mmax(1);
				fs_new = [fs_new fs(k)+mmax];
			end
		end
		fs_new = [fs_new fs(end)];
	end
	Foptj = Fn(fs(2),fs(1));
	for k=2:length(fs)-1
		Foptj = Foptj * Fn(fs(k+1),fs(k));
	end
	pot.Fopt(j) = Foptj; % fitness of optimal sequence of generations
	pot.tgen{j} = t1(fs(1:end-1)); % start times of each generation
	pot.ngen(j) = (length(fs)-1)/p.Nyears; % generations per year
end
v.pot = pot;

% save a few scalar stats at the level where coltraneEnsemble.m will find them
v.Fpot = max(v.pot.Fopt); % total eggs per egg over the model run
f = find(v.pot.Fopt==max(v.pot.Fopt));
f = round(mean(f));
v.tFpot = p.tspawn(f);
v.Sapot= max(v.S(:,f));
v.ngenpot = v.pot.ngen(f); % generations per year for optimal strategy
v.Fpot_per_yr = v.Fpot ^ (1/p.Nyears); % eggs per egg per year
v.Fpot_per_gen = v.Fpot ^ (1/p.Nyears/v.ngenpot); % eggs per egg per generation
	% (averaged over generations in the optimal strategy)
v.tgen1pot = min(yearday([v.pot.tgen{f}]));
	% first spawning yearday in the optimal strategy
v.afracpot = sum(v.a(:,f)==1)/sum(~isnan(v.a(:,f)));
	% amount of time spent in diapause in the optimal strategy
	
% stats for strictly 1 and 2 yrs/gen strategies
v.Fpot_at1yr = max(v.pot.Fmult(1,:)); 
v.tFpot_at1yr = nan;
v.Sapot_at1yr = nan;
if ~isnan(v.Fpot_at1yr)
	f = find(v.pot.Fmult(1,:)==max(v.pot.Fmult(1,:)));
	f = round(mean(f));
	v.tFpot_at1yr = p.tspawn(f);
	v.Sapot_at1yr= max(v.S(:,f));
end	
v.Fpot_at2yr = nan;
v.tFpot_at2yr = nan;
v.Sapot_at2yr = nan;
if p.Nyears > 1
	v.Fpot_at2yr = max(v.pot.Fmult(2,:));
	if ~isnan(v.Fpot_at2yr)
		f = find(v.pot.Fmult(2,:)==max(v.pot.Fmult(2,:)));
		f = round(mean(f));
		v.tFpot_at2yr = p.tspawn(f);
		v.Sapot_at2yr= max(v.S(:,f));
	end
end


% ------------------------------------------------------------------------------
% egg/reserve model postprocessing ---------------------------------------------

% metrics for each cohort
% (v.by_t.* are raw metrics for each spawning date; v.* are metrics for the
% population, once a stable annual cycle of spawning has been determined)

% actual adult size
Na = exp(v.lnN);
SR = v.S+v.R;
Na(v.D<1 | isnan(SR)) = 0;
SR(isnan(SR)) = 0;
v.by_t.Wa = sum(SR.*Na)./sum(Na); % abundance-weighted average of adult states
v.by_t.We = p.r_ea .* v.by_t.Wa.^p.exp_ea;
	% updated egg size: if Wa = p.Wa0, We = p.We0

% lifetime egg production (1-generation calculation)
v.LEP = sum((v.Einc+v.Ecap).*exp(v.lnN))./v.by_t.We.*p.dt;
v.LEP(isnan(v.LEP)) = 0;
% note: a simple way to take into account internal life-history mismatch without
% the eigenvalue calculation below would be a two-generation calculation:
% instead of LEP = integral of E * N, LEP2 = integral of E * N * LEP

% terms in the lifetime energy budget
% yield to predators = total gain - egg prod beyond replacement
de = (v.Fin-v.Fmet).*v.S.*exp(v.lnN);
de(isnan(de)) = 0;
v.by_t.ener_gain = p.dt .* sum(de);
v.by_t.ener_egg = (v.LEP - 1).*v.by_t.We;
v.by_t.ener_yield = v.by_t.ener_gain - v.by_t.ener_egg; 
de = v.m.*(v.S+v.R).*exp(v.lnN);
de(isnan(de)) = 0;
v.by_t.ener_yield1 = p.dt .* sum(de); % for checking
% there is in fact a difference between yield = gain - eggs and yield1 = 
% explicit integral of mortality, which comes from the final value
% of (R+S)N at the moment of starvation: net accumulation of biomass in the
% cohort itself. But this might well be thought of as (eventually) another
% yield to predators.
de = v.m.*v.R.*exp(v.lnN);
de(isnan(de)) = 0;
v.by_t.ener_yield1_R =  p.dt .* sum(de);
v.by_t.rhoY = v.by_t.ener_yield1_R ./ v.by_t.ener_yield1;
	% yield (trophic transfer) of R as fraction of total yield

% other metrics by cohort
v.by_t.develtime = sum(v.D < 1) .* p.dt;
D1 = v.D;
D1(v.a~=0) = nan;
v.by_t.Ddiamin = min(D1);
v.rho = max(0,v.R) ./ (v.R + v.S);
WN = (v.R+v.S).*exp(v.lnN);
WN(isnan(WN))=0;
rhoWN = v.rho .* WN;
rhoWN(isnan(rhoWN))=0;
v.by_t.rhomean = sum(rhoWN) ./ sum(WN);
Pgrow = forcing.P;
Pgrow(v.a < 1 | v.D >=1) = nan;
v.by_t.satmean = nanmean(Pgrow ./ (p.Ks + Pgrow));
v.by_t.capfrac = nanmean(v.Ecap.*exp(v.lnN)) ./ ...
				 nanmean((v.Ecap+v.Einc).*exp(v.lnN));
NA = exp(v.lnN);
NA(v.D<1) = 0;
v.by_t.recruitmentToAdult = max(NA);


% find a stable cycle of egg production

% V0 ~ EN/We dt; sum over first dimension (time) is LEP
v.V0 = (v.Einc+v.Ecap).*exp(v.lnN)./repmat(v.by_t.We,[N 1]).*p.dt;
v.V0(isnan(v.V0)) = 0;
% chunk it into an annual cycle (by tspawn, ignoring which year a given bit of 
% spawning happens). Make sure p.dt is a divisor of p.dt_tspawn, and both fit
% evenly into a 365-day year!
v.V = zeros(NC,NC);
i0 = repmat((1:(p.dt_tspawn/p.dt))',[1 p.Nyears]);
i0 = i0 + repmat((NC-1).*size(i0,1).*(0:p.Nyears-1),[size(i0,1) 1]);
i0 = i0(:);
for j=1:NC
	v.V(j,:) = v.V(j,:) + sum(v.V0(i0 + (j-1)*(p.dt_tspawn/p.dt), :));
end
v.by_t.genlength = sum((forcing.t - repmat(p.tspawn(:)',[N 1])) .* v.V0) ...
				./ sum(v.V0); % time from tspawn to mean egg prod date

% v.V is a transition matrix for generating an annual time series of eggs 
% spawned from the previous year of eggs spawned. The first eigenvector of v.V 
% is thus a stable annual cycle of egg production.
% n(t)_next = V' * n(t)'
[evecs,evals] = eig(v.V');
eval1 = max(real(diag(evals)));
if eval1 > 0
	v.lam = eval1 - 1; % first eigenvalue is 1 + population growth rate lambda
	f = find(real(diag(evals))==eval1);
	v.n = evecs(:,f(1))';
	v.n(isnan(v.n)) = 0;
	v.n = v.n ./ sum(v.n);
		% normalize so that v.n is _relative_ distribution of spawning dates--
		% a set of weights, and therefore proportional to p.dt_tspawn
else
	v.lam = nan;
	v.n = zeros(size(p.tspawn));
end
	
% population-level metrics and cycles
fields = fieldnames(v.by_t);
for i=1:length(fields)
	F = v.by_t.(fields{i});
	F(isnan(F))=0;
	v.(fields{i}) = sum(F.*v.n.*v.LEP) ./ sum(v.n.*v.LEP);
		% weighted by egg production and egg fitness (expectation value in next 
		% gen). To weight these simply by egg production, sum(F.*v.n)
end
tt = [p.tspawn(cumsum(v.n)>0.1.*sum(v.n)) nan];
v.tE10 = tt(1); % yearday when 1st 10% of egg prod completed
tt = [p.tspawn(cumsum(v.n)>0.5.*sum(v.n)) nan];
v.tE50 = tt(1);
tt = [p.tspawn(cumsum(v.n)>0.9.*sum(v.n)) nan];
v.tE90 = tt(1);
v.tEmax = nanmean(p.tspawn(v.n==max(v.n)));
v.LEPn = sum(v.LEP .* v.n); % egg-production-weighted LEP (compare with lam)


% abundance-weights for adding up time series to make a population.
% Nn = relative spawning intensity n * survivorship of each case N, taking into 
% account carryover from previous years (this last point is what makes the code 
% a bit complicated)
kk = 1:(p.dt_tspawn/p.dt):(365/p.dt);
rsh = [365/p.dt p.Nyears+1];
Nn0 = exp(v.lnN) .* repmat(v.n,[N 1]);
Nntot = sum(sum(reshape(Nn0(1:end-1,:)',[NC rsh]),3),1);
Nn = [Nn0(1:end-1,:) ./ repmat(Nntot(:),[p.Nyears+1 NC]); zeros(1,NC)];
% another version, for diapause-capable stages only
Nndia0 = exp(v.lnN) .* repmat(v.n,[N 1]);
Nndia0(v.D < p.Ddia) = 0;
Nndiatot = sum(sum(reshape(Nndia0(1:end-1,:)',[NC rsh]),3),1);
Nndia = [Nndia0(1:end-1,:) ./ repmat(Nndiatot(:),[p.Nyears+1 NC]); zeros(1,NC)];
% another version, biomass-weighted
WNn0 = (v.R+v.S) .* Nn0;
WNn0(isnan(WNn0)) = 0;
WNntot = sum(sum(reshape(WNn0(1:end-1,:)',[NC rsh]),3),1);
WNn = [WNn0(1:end-1,:) ./ repmat(WNntot(:),[p.Nyears+1 NC]); zeros(1,NC)];
% stage structure, abundance-weighted (_N) and biomass-weighted (_W)
stages = 1:13;
stage = D2stage(v.D);
stageFracN0 = nan.*ones(N,length(stages));
stageFracW0 = nan.*ones(N,length(stages));
for i=1:length(stages)
	stageFracN0(:,i) = sum((stage==stages(i)) .* Nn, 2);
	stageFracW0(:,i) = sum((stage==stages(i)) .* WNn, 2);
end
stageFracN0 = sum(reshape(stageFracN0(1:end-1,:)',[length(stages) rsh]),3)';
stageFracW0 = sum(reshape(stageFracW0(1:end-1,:)',[length(stages) rsh]),3)';
v.stageFracN = stageFracN0(kk,:);
v.stageFracW = stageFracW0(kk,:);
% other population time series
afrac0 = sum((v.a==1) .* Nndia, 2);
afrac0 = sum(reshape(afrac0(1:end-1),rsh),2);
v.afrac = afrac0(kk);
	% fraction of diapause-capable individuals that are in diapause
