function OUT = coltrane_integrate(forcing,p,t0,s,whatToSave);

%   v = coltrane_integrate(forcing,p,t0,s,'everything');
% dF1 = coltrane_integrate(forcing,p,t0,s,'fitness only');
%
% calculates a(t), D(t), W(t), N(t), E(t), and dF1 for set of spawning dates t0
% and a set of strategies s. The fields of s should be tdia_exit, tdia_enter,
% and dtegg.
%
% the user doesn't call this; coltraneModel.m does.


if nargin < 5, whatToSave = 'everything'; end

v = forcing; % copy over forcing time series
NT = size(v.t,1); % # timesteps
% determine the simulation timestep from the forcing time series
dt = v.t(2) - v.t(1);

% get t, t0, and s to the proper shapes
NC = length(t0);
fields = fieldnames(s);
NS = length(s.(fields{1})(:));
v.t = repmat(v.t(:), [1 NC NS]);
v.t0 = repmat(reshape(t0,[1 NC 1]), [1 1 NS]);
for k=1:length(fields)
	v.(fields{k}) = repmat(reshape(s.(fields{k}),[1 1 NS]),[1 NC 1]);
end
tegg = v.t0 + v.dtegg; % date that egg production starts, [1 NC NS]
tegg_3d = repmat(tegg,[NT 1 1]);

dF1 = zeros(NT,NC,NS);

% activity: a(t) ------------------------------------------------------
v.a = double(~(v.yday >= repmat(v.tdia_enter,[NT 1 1]) | ...
			   v.yday <= repmat(v.tdia_exit,[NT 1 1])));
	% assumes that tdia_enter > tdia_exit
t0yday = reshape(yearday(v.t0),size(v.t0));
activeSpawning = ~(t0yday >= v.tdia_enter | t0yday <= v.tdia_exit);
								
% temperature response factors: qd, qg  --------------------------------
v.temp = v.a .* v.T0 + (1-v.a) .* v.Td;
qd = (p.Q10d^0.1) .^ v.temp;
qg = (p.Q10g^0.1) .^ v.temp;

% prey saturation ------------------------------------------------------
v = preySaturation(v,p);
if size(v.sat,2)==1
	v.sat = repmat(v.sat(:),[1 NC NS]);
end

% development: D(t) ----------------------------------------------------
isalive = v.t >= repmat(v.t0,[NT 1 1]); % spawned yet?
if p.requireActiveSpawning
	% eliminate (t0,s) combinations in which t0 falls during diapause
	isalive = isalive & repmat(activeSpawning,[NT 1 1]); 
end
dDdt = isalive .* p.u0 .* qd; % nonfeeding formula
v.D = cumsum(dDdt) .* dt;
isfeeding = v.D >= p.Df;
dDdt_feeding = isalive .* p.u0 .* qd .* v.a .* v.sat;
dDdt(isfeeding) = dDdt_feeding(isfeeding); % full formula 
v.D = cumsum(dDdt) .* dt;
v.D(v.D>1) = 1;

% the life history is divided into two phases, growth and egg production. 
% This is equivalent to juvenile and adult if the animals delay their
% final maturation until just before egg prod. begins, but not if they enter
% C6 and then delay egg prod. (e.g. in the case of overwintering C6's)
isineggprod = v.t >= tegg_3d;
isgrowing = v.D >= p.Df & ~isineggprod;

% check D(t) for consistency with dtegg
toolate = any(v.D<1 & isineggprod);
if all(toolate)
	v.D = [];
	if strcmpi(whatToSave,'fitness only'), OUT = dF1; else, OUT = v; end
	return
	% if there are no good solutions at all, don't bother with the rest of the
	% model. isempty(v.D) is the test of whether the model completed
else
	v.D(:,toolate) = nan;
	% if there's an inconsistency for some (t0,tdia_enter,tdia_exit) and not
	% others, blank out the bad ones 
end

% calculate D in middle of first winter, just as a diagnostic
[yr0,~] = datevec(v.t0);
first31dec = repmat(reshape(datenum(yr0,12,31), [1 NC NS]), [NT 1 1]);
is31dec = abs(v.t-first31dec) == ...
			  repmat(min(abs(v.t-first31dec)),[NT 1 1]);
v.D_winter = reshape(v.D(is31dec),[1 NC NS]);

% flag time points at which the animal is in diapause but at a 
% diapause-incapable stage (D < Ds in this version: Ddia has been dropped)
isactive = isalive & v.a==1;
hasbeenactive = (cumsum(isactive) >= 1);
isfailingtodiapause = isalive & hasbeenactive & v.a==0 & v.D < p.Ddia;
hasfailedtodiapause = (cumsum(isfailingtodiapause)>1);
v.D(hasfailedtodiapause) = nan;

isalive = isalive & ~isnan(v.D);
isineggprod = isineggprod & isalive;
isgrowing = isgrowing & isalive;

% pick a guess at adult size based on mean temperature in the forcing,
% and pick a corresponding guess at egg size based on this.
% in v1.0 these were called p.Wa0 and p.We0.
T_nominal = mean(v.temp);
qoverq = (p.Q10g./p.Q10d).^(T_nominal./10);
co = (1-p.theta) .* p.GGE_nominal .* (1-p.Df) .* qoverq;
v.Wa_theo = (co .* p.I0 ./ p.u0) .^ (1./(1-p.theta));
v.We_theo = p.r_ea .* v.Wa_theo .^ p.exp_ea;

% energy gain and growth: G(t), W(t) -------------------------------------------
% approximate curve of W(D), valid for a=1
% (this is used for the allometry in the calculation of net gain)
satmean = sum(v.sat .* isgrowing) ./ sum(isgrowing);
c = (1-p.theta) .* qg ./ qd .* p.I0 ./ p.u0 .* ...
	(p.r_assim - p.rm ./ repmat(satmean,[NT 1 1 1]));
c = max(c,0);
Wapprox = zeros(size(v.D));
We_theo_3d = repmat(v.We_theo,[NT 1 1]);
Wapprox(isgrowing) = ((We_theo_3d(isgrowing).^(1-p.theta)) + ...
	c(isgrowing) .* (v.D(isgrowing) - p.Df)) .^ (1/(1-p.theta));
% net gain during the growth phase (here stored as gain * W)
ImaxW = zeros(size(Wapprox));
ImaxW(isgrowing) = ...
	qg(isgrowing) .* p.I0 .* Wapprox(isgrowing) .^ p.theta;
astar = p.rb + (1-p.rb) .* v.a;
GW = zeros(size(v.D));
GW(isgrowing) = ImaxW(isgrowing) .* ...
	(v.a(isgrowing) .* p.r_assim .* v.sat(isgrowing) - ...
	p.rm .* astar(isgrowing));
% integrate to get the correct W over development
v.W = We_theo_3d + cumsum(GW) .* dt;
v.W = max(0,v.W);
% calculate the fraction of this which is reserves
fs = ones(size(v.W));
isstoringR = v.D >= p.Ds;
fs(isstoringR) = (1-v.D(isstoringR)) ./ (1-p.Ds);
fs(GW < 0) = 0;
v.R = cumsum((1-fs).*GW) .* dt;
% adult size Wa (= size at the moment egg prod begins)
last = ~isineggprod(1:end-1,:,:) & isineggprod(2:end,:,:);
v.Wa = max(v.W(1:end-1,:,:) .* last);
Wa_3d = repmat(v.Wa,[NT 1 1]);
v.W(isineggprod) = Wa_3d(isineggprod);
	% at this point, set W(t) = Wa throughout the egg-production phase;
	% once we calculate Ecap(t) below, we will subtract it
% likewise for reserves
v.Ra = max(v.R(1:end-1,:,:) .* last); % the max() shouldn't be necessary
Ra_3d = repmat(v.Ra,[NT 1 1]);
v.R(isineggprod) = Ra_3d(isineggprod);
% calculate net gain for the egg-prod phase
% (and for the growth phase as well: previously we calculated GW but not G)
Imax = qg .* p.I0 .* v.W.^(p.theta-1);
Imax(~isfeeding) = 0;
Imax(v.W < We_theo_3d) = 0;
I = v.a .* p.r_assim .* v.sat .* Imax;
	% this was r_assim * I, as opposed to I, in Coltrane 1.0
M = p.rm .* astar .* Imax;
v.G = I - M;

% check for starvation
isstarving = isgrowing & v.R < -p.rstarv .* v.W;
isalive = isalive & cumsum(isstarving)==0;
isineggprod = isineggprod & isalive;

% mortality and survivorship: N(t) -------------------------------------
v.m = p.m0 .* qg .* v.a .* v.W.^(p.theta-1); % mort. rate at T, size
v.m(~isalive) = 0;
v.lnN = cumsum(-v.m) .* dt;
% calculate adult recruitment (= recruitment at the moment egg prod begins)
lnNa = v.lnN;
lnNa(~isineggprod) = nan;
v.Na = exp(max(lnNa));

% egg production -------------------------------------------------------
% income egg prod = GW once tegg is reached
v.Einc = max(0, v.G .* v.W .* isineggprod);
if p.requireActiveSpawning
	v.Einc = v.Einc .* v.a;
end
% start by calculating what capital egg prod would be if there were
% inexhaustible energy for it
Emax = zeros(size(v.G));
Emax(isineggprod) = Imax(isineggprod) .* v.W(isineggprod);
if p.requireActiveSpawning
	Emax = Emax .* v.a;
end
Ecap = max(0, Emax - v.Einc);
% debit this from reserves until the point where it exceeds what's available
% (may have an error up to Imax * dt)
deltaR = cumsum(Ecap).*dt;
toomuch = deltaR > Ra_3d;
Ecap(toomuch) = 0;
deltaR = cumsum(Ecap).*dt;
v.R = v.R - deltaR;
v.W = v.W - deltaR;
% add them up
v.E = (v.Einc + Ecap);
v.capfrac = sum(Ecap) ./ sum(v.E);

% contributions to fitness at each t -------------------------------------------
v.dF1 = (v.E .* exp(v.lnN) .* dt) ./ v.We_theo;
v.dF1(isnan(v.dF1)) = 0;
		
% clean up -------------------------------------------------------------
if strcmpi(whatToSave,'fitness only')
	OUT = v.dF1;
else
	v.a(~isalive) = nan;
	v.D(~isalive) = nan;
	v.W(~isalive) = nan;
	v.R(~isalive) = nan;
	v.lnN(~isalive) = -Inf;
	v.E(~isalive) = 0;
	OUT = v;
end
