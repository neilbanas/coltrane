function OUT = coltrane_integrate(forcing,p,t0,s,whatToSave)

%   v = coltrane_integrate(forcing,p,t0,s,'everything');
%   								  ...,'scalars only');
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
tayday = reshape(yearday(v.t0+v.dtegg),size(v.t0));
activeSpawning2 = ~(tayday >= v.tdia_enter | tayday <= v.tdia_exit);
								
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
isineggprod = v.t >= tegg_3d;
if p.requireActiveSpawning
	% eliminate (t0,s) combinations in which t0 falls during diapause
	isalive = isalive & repmat(activeSpawning,[NT 1 1]); 
	% likewise (t0,s) combinations in which t0+dtegg falls during diapause
	isalive = isalive & repmat(activeSpawning2,[NT 1 1]);
    isineggprod = isineggprod & v.a;
end
dDdt = isalive .* p.u0 .* qd; % nonfeeding formula
v.D = cumsum(dDdt) .* dt;
isfeeding = v.D >= p.Df;
dDdt_feeding = isalive .* p.u0 .* qd .* v.a .* v.sat;
dDdt(isfeeding) = dDdt_feeding(isfeeding); % full formula 
v.D = cumsum(dDdt) .* dt;
v.D(v.D>1) = 1;
isfeeding = v.D >= p.Df; % recalculated because cumsum() is one timestep
						 % off from a true forward integration
                         
% the life history is divided into two phases, growth and egg production. 
% This is equivalent to juvenile and adult if the animals delay their
% final maturation until just before egg prod. begins, but not if they enter
% C6 and then delay egg prod. (e.g. in the case of overwintering C6's)
% isineggprod = v.t >= tegg_3d;
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
% diapause-incapable stage, and mark these cases as dead
isactive = isalive & v.a==1;
hasbeenactive = (cumsum(isactive) >= 1);
isfailingtodiapause = isalive & hasbeenactive & v.a==0 & v.D < p.Ddia;
hasfailedtodiapause = (cumsum(isfailingtodiapause)>1);
v.D(hasfailedtodiapause) = nan;

isalive = isalive & ~isnan(v.D);
isfeeding = isfeeding & isalive;
isgrowing = isgrowing & isalive;
isineggprod = isineggprod & isalive;
fullydeveloped = v.D == 1; % distinct from isineggprod

% pick a guess at adult size based on mean temperature in the forcing,
% and pick a corresponding guess at egg size based on this.
% in v1.0 these were called p.Wa0 and p.We0.
T_nominal = mean(v.temp);
qoverq = (p.Q10g./p.Q10d).^(T_nominal./10);
co = (1-p.theta) .* p.GGE_nominal .* (1-p.Df) .* qoverq;
v.Wa_theo = (co .* p.I0 ./ p.u0) .^ (1./(1-p.theta));
v.We_theo = p.r_ea .* v.Wa_theo .^ p.exp_ea;

% energy gain, growth, egg production: G(t), W(t), R(t), E(t) ------------------
v.G = zeros(size(v.D));
v.W = repmat(v.We_theo,[NT 1 1]);
v.R = zeros(size(v.D));
v.Einc = zeros(size(v.D));
v.E = zeros(size(v.D));
astar = p.rb + (1-p.rb) .* v.a;

% track maximum reserves to limit growth at C6
if NS > 1 || NS == 0
    maxReserves = zeros(NC,NS);
else
    maxReserves = zeros(1,NC);
end

for n=1:NT-1
	f = isfeeding(n,:); g = isgrowing(n,:);
    d = fullydeveloped(n,:); e = isineggprod(n,:);
	% net gain
	Imax_nf = qg(n,f) .* p.I0 .* v.W(n,f).^(p.theta-1);
	I_nf = v.a(n,f) .* p.r_assim .* v.sat(n,f) .* Imax_nf;
	M_nf = p.rm .* astar(n,f) .* Imax_nf;
	v.G(n,f) = I_nf - M_nf;
    GWdt = v.G(n,:) .* v.W(n,:) .* dt;
    energyLost = GWdt < 0;
    gdGain = g & d & ~energyLost; % growth phase, D=1, energy gained
    eGain = e & ~energyLost; % egg phase, energy gained
    % allocation to growth
	dW = GWdt;
    dW(gdGain) = min(dW(gdGain), max(0, maxReserves(gdGain) - v.R(n,gdGain))); % restricted weight gain when fully developed in growth phase
    dW(eGain) = 0; % egg phase energy gains do not increase weight    
    v.W(n+1,:) = max(0, v.W(n,:) + dW);
	% allocation to reserves
	fr = (v.D(n,:) - p.Ds) ./ (1 - p.Ds);
	fr = max(0,min(1,fr));
	fr(energyLost) = 1; % all net losses come from R
    v.R(n+1,:) = v.R(n,:) + fr .* dW;
    gainR = squeeze(v.R(n+1,:,:)) > maxReserves;
    maxReserves(gainR) = v.R(n+1,gainR);    
    % income egg production
    v.Einc(n,:) = max(0,GWdt)./dt .* e;
	% capital egg production
	Emax = zeros(size(GWdt));
	Emax(f) = p.r_assim .* Imax_nf .* e(f) .* v.W(n+1,f);
    Ecap = max(0, Emax - v.Einc(n,:));
    dR = min(max(0,v.R(n+1,:)), Ecap.*dt);
    v.E(n,:) = v.Einc(n,:) + dR./dt;
	v.W(n+1,:) = v.W(n+1,:) - dR;
	v.R(n+1,:) = v.R(n+1,:) - dR;
end


% adult size Wa, Ra (= size at the moment egg prod begins)
last = ~isineggprod(1:end-1,:,:) & isineggprod(2:end,:,:);
v.Wa = max(v.W(1:end-1,:,:) .* last);
v.Ra = max(v.R(1:end-1,:,:) .* last); % the max() shouldn't be necessary
% update the estimate of We
v.We = p.r_ea .* v.Wa .^ p.exp_ea;
% capital fraction of egg production
v.capfrac = sum(v.E - v.Einc) ./ sum(v.E);

% check for starvation
isstarving = v.R < -p.rstarv .* v.W;
isalive = isalive & cumsum(isstarving)==0;
isineggprod = isineggprod & isalive;
v.E(~isineggprod) = 0;

% mortality and survivorship: N(t) ---------------------------------------------
v.m = p.m0 .* qg .* v.a .* v.W.^(p.theta-1); % mort. rate at T, size
v.m(~isalive) = 0;
v.lnN = cumsum(-v.m) .* dt;
% calculate adult recruitment (= recruitment at the moment egg prod begins)
lnNa = v.lnN;
lnNa(~isineggprod) = nan;
v.Na = exp(max(lnNa));

% contributions to fitness at each t -------------------------------------------
% v.dF1 = real(v.E .* exp(v.lnN) .* dt) ./ v.We_theo;
v.dF1 = (v.E .* exp(v.lnN)) ./ v.We_theo .* dt;
v.dF1(isnan(v.dF1)) = 0;
v.tEcen = sum(v.t .* v.dF1) ./ sum(v.dF1); % center of mass of E*N
	% this minus t0 is generation length
    
% days between reaching maximum possible size and
% t0+dtegg, when energy gains are wasted as they cant go to reserve gorowth
% or income egg production.
% maxSize = v.R ./ v.W >= p.maxRelReserves;
% v.noGrowthDaysPreDia = sum(maxSize & cumsum(~isactive) == 0);
% v.noGrowthDaysPostDia = sum(maxSize & cumsum(isactive) == 2);

% clean up ---------------------------------------------------------------------
if strcmpi(whatToSave,'fitness only')
	OUT = v.dF1;
elseif strcmpi(whatToSave,'scalars only')
	OUT = keepScalars(v);
else
	v.a(~isalive) = nan;
	v.D(~isalive) = nan;
	v.W(~isalive) = nan;
	v.R(~isalive) = nan;
	v.lnN(~isalive) = -Inf;
	v.E(~isalive) = 0;
	OUT = v;
end
