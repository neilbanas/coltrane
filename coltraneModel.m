function OUT = coltraneModel(forcing0,p);

% out = coltraneModel(forcing0);
% ... = coltraneModel(forcing0,p);
%
% time-integrates one population with a single set of traits, but variable
% spawning date tspawn. Calculates summary metrics for both the phi and ER
% versions of the model.
%
% For full model description, see Banas et al (2016) Front. Mar. Res., submitted
% (http://neilbanas.com/projects/coltrane)


if nargin<2, p = coltraneParams; end
if nargin<1, forcing0 = coltraneForcing('simple-g'); end

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
v.S = nan.*ones(N,NC); % structural biomass
v.phi = nan.*ones(N,NC); % potential reserves/eggs
v.a = nan.*ones(N,NC); % activity level
v.R = nan.*ones(N,NC); % reserve biomass
v.Einc = zeros(N,NC); % income egg production
v.Ecap = zeros(N,NC); % capital egg production
v.m = zeros(N,NC); % mortality rate
v.Fin = zeros(N,NC); % ingestion flux
v.Fmet = zeros(N,NC); % metabolism flux

% pick a guess at adult size (Wa0) based on mean temperature in the forcing,
% and pick a corresponding guess at egg size (We0) based on this.
% these are no longer calculated in coltraneParams.m but still stored in _p_
% until a bigger reorganisation to come
T_nominal = mean(forcing.T(:,1));
qoverq = (p.Q10g./p.Q10d).^(T_nominal./10);
co = (1-p.theta) .* p.GGE_nominal .* (1-p.Df) .* qoverq;
p.Wa0 = (co .* p.I0 ./ p.u0) .^ (1./(1-p.theta));
p.We0 = p.r_ea .* p.Wa0.^p.exp_ea;

% initial conditions
for i=1:NC
	f = find(t == v.tspawn(i));	
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
	sp = forcing.t(n,:) >= v.tspawn; % the ones that have been spawned so far
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
	fs = (min(v.D(n,:),1)-1) ./ (p.Ds - 1 - small);
	fs = min(1, max(0, fs)); % fs=1 means entirely to/from S
	GSdtpos = max(0,GSdt);
	v.S(n+1,sp) = v.S(n,sp) + fs(sp) .* GSdtpos(sp);
	v.phi(n+1,sp) = v.phi(n,sp) + (1-fs(sp)) .* GSdtpos(sp);
	% negative net gain is taken from phi
	GSdtneg = min(0,GSdt);
	v.phi(n+1,sp) = v.phi(n+1,sp) + GSdtneg(sp);	
	% determine diapause
	if p.myopicDiapause
		% instantaneous assessment; fluctuates step by step
		C = max(0, 1 + min(p.r_phi_max, v.phi(n,:))./v.S(n,:));
		sat_crit = p.rm .* (1-p.rb) ./ p.r_assim + C .* p.m0./p.I0./p.r_assim;
		dia = v.D(n,:) > p.Ddia & sat < sat_crit;
	else
		% an explicit range of yeardays
		yday = mod(forcing.t(n,:),365);
		if p.tdia_enter > p.tdia_exit
			dia = v.D(n,:) > p.Ddia & ...
				(yday >= p.tdia_enter | yday <= p.tdia_exit);
		else
			dia = v.D(n,:) > p.Ddia & ...
				(yday >= p.tdia_enter & yday <= p.tdia_exit);
		end 
	end
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
	v.R(n+1,sp) = v.R(n+1,sp) - v.Einc(n,sp).*p.dt;
	v.Ecap(n,sp & egg) = p.fcap .* min(max(0,v.R(n+1,sp & egg)./p.dt), ...
						           Emax(sp & egg) - v.Einc(n,sp & egg));
	v.R(n+1,sp) = v.R(n+1,sp) - v.Ecap(n,sp).*p.dt;
	
	% save fluxes for use in energy budget
	v.Fin(n,:) = Fin;
	v.Fmet(n,:) = Fmet;
end


v.t = forcing.t(:,1);
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


% organise output into structures ----------------------------------------------
OUT.p = p;
OUT.forcing = forcing;
OUT.var = v;


% postprocessing ---------------------------------------------------------------
OUT.pot = coltranePostproc_phi(v,p,forcing);
[OUT.cohort, OUT.pop, OUT.routine] = coltranePostproc_ER(v,p,forcing);

