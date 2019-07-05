function v = coltrane_integrate(forcing,p,t0,s,whatToSave);

% v = coltrane_integrate(forcing,p,t0,s,'everything');
%   							    ...,'scalars only');
%   								...,'scalars and fitness');
%
% calculates a(t), D(t), W(t), N(t), E(t), and dF1 for a set of spawning dates
% t0 and a set of strategies s.
% 
% The fields of s should be tdia_exit, tdia_enter, and dtegg. They can be any
% shape; they will be rearranged into a flat list.


if nargin < 5, whatToSave = 'everything'; end

v = forcing; % copy over forcing time series
NT = size(v.t,1); % # timesteps
v.yday = yearday(v.t); % make sure this is filled in and consistent
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

v.F1 = zeros(1,NC,NS);
v.F2 = zeros(1,NC,NS);

% activity: a(t) ------------------------------------------------------
v.a = double(~(v.yday >= repmat(v.tdia_enter,[NT 1 1]) | ...
			   v.yday <= repmat(v.tdia_exit,[NT 1 1])));
	% assumes that tdia_enter > tdia_exit
isalive = v.t >= repmat(v.t0,[NT 1 1]); % born yet?
if p.requireActiveSpawning
	% eliminate (t0,s) combinations in which t0 falls during diapause
	t0_yday = yearday(v.t0);
	activeSpawning = ~(t0_yday >= v.tdia_enter | t0_yday <= v.tdia_exit);
	isalive = isalive & repmat(activeSpawning,[NT 1 1]); 
	% likewise (t0,s) combinations in which t0+dtegg falls during diapause
	% (this eliminates some redundancy)
	ta_yday = yearday(v.t0+v.dtegg);
	activeNextSpawning = ~(ta_yday >= v.tdia_enter | ta_yday <= v.tdia_exit);
	isalive = isalive & repmat(activeNextSpawning,[NT 1 1]); 
end


% temperature response factors: qd, qg  --------------------------------
v.temp = v.a .* v.T0 + (1-v.a) .* v.Td;
qd = zeros(NT,NC,NS);
qd(isalive) = (p.Q10d^0.1) .^ v.temp(isalive);
qg = zeros(NT,NC,NS);
qg(isalive) = (p.Q10g^0.1) .^ v.temp(isalive);


% prey saturation ------------------------------------------------------
v = preySaturation(v,p);
if size(v.sat,2)==1
	v.sat = repmat(v.sat(:),[1 NC NS]);
end


% development: D(t) ----------------------------------------------------
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
isineggprod = v.t >= tegg_3d;

% check D(t) for consistency with dtegg
toolate = any(v.D<1 & isineggprod);
if all(toolate)
	v.D = [];
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
is31dec = is31dec & cumsum(is31dec)==1;
v.D_winter = reshape(v.D(is31dec),[1 NC NS]);

% flag time points at which the animal is in diapause but at a 
% diapause-incapable stage, and mark these cases as dead
isactive = isalive & (v.a==1 | ~isfeeding);
hasbeenactive = (cumsum(isactive) >= 1);
isfailingtodiapause = isalive & hasbeenactive & v.a==0 & v.D < p.Ddia;
hasfailedtodiapause = (cumsum(isfailingtodiapause)>1);
v.D(hasfailedtodiapause) = nan;

isalive = isalive & ~isnan(v.D);
isfeeding = isfeeding & isalive;
isineggprod = isineggprod & isalive;


% pick a guess at adult size based on mean temperature in the forcing,
% and pick a corresponding guess at egg size based on this.
% in v1.0 these were called p.Wa0 and p.We0.
T_nominal = mean(v.temp);
qoverq = (p.Q10g./p.Q10d).^(T_nominal./10);
co = (1-p.theta) .* p.GGE_nominal .* (1-p.Df) .* qoverq;
v.Wa_theo = (co .* p.I0 ./ p.u0) .^ (1./(1-p.theta));
v.We_theo = p.r_ea .* v.Wa_theo .^ p.exp_ea;


% energy gain, growth, egg production: G(t), W(t), R(t), E(t) ------------------
v.G = zeros(NT,NC,NS);
v.W = repmat(v.We_theo,[NT 1 1]);
v.R = zeros(NT,NC,NS);
v.Einc = zeros(NT,NC,NS);
v.E = zeros(NT,NC,NS);
astar = p.rb + (1-p.rb) .* v.a;

for n=1:NT-1
	f = isfeeding(n,:); e = isineggprod(n,:);
	% net gain
	Imax_nf = qg(n,f) .* p.I0 .* v.W(n,f).^(p.theta-1);
	I_nf = v.a(n,f) .* p.r_assim .* v.sat(n,f) .* Imax_nf;
	M_nf = p.rm .* astar(n,f) .* Imax_nf;
	v.G(n,f) = I_nf - M_nf;
    GWdt = v.G(n,:) .* v.W(n,:) .* dt;
	% allocation to growth
	dW = GWdt;
	dW(dW>0 & e) = 0; % definitely
%	dW(dW>0 & v.D(n,:)==1) = 0; % maybe
		% this is meant to accomplish what Aidan's maxReserves mechanism did.
		% note that positive gain is simply lost between D=1 and the start of egg prod,
		% as if the animals are wishing they were in diapause. It might be too strict.	
	v.W(n+1,:) = max(0, v.W(n,:) + dW);
	% allocation to reserves
	fr = (v.D(n,:) - p.Ds) ./ (1 - p.Ds);
	fr = max(0,min(1,fr));
	fr(GWdt < 0) = 1; % all net losses come from R
    v.R(n+1,:) = v.R(n,:) + fr .* dW;
    % income egg production
    v.Einc(n,:) = max(0,GWdt)./dt .* e;
	% capital egg production
	Emax = zeros(size(GWdt));
	Emax(f) = p.r_assim .* Imax_nf .* e(f) .* v.W(n,f);
    Ecap = max(0, Emax - v.Einc(n,:));
    dR = min(max(0,v.R(n+1,:)), Ecap.*dt);
    v.E(n,:) = v.Einc(n,:) + dR./dt;
	v.W(n+1,:) = v.W(n+1,:) - dR;
	v.R(n+1,:) = v.R(n+1,:) - dR;
end

% adult size Wa, Ra (= size at the moment egg prod begins)
last = ~isineggprod(1:end-1,:,:) & isineggprod(2:end,:,:);
v.Wa = max(v.W(1:end-1,:,:) .* last);
v.Ra = max(v.R(1:end-1,:,:) .* last);
% update the estimate of We
v.We = p.r_ea .* v.Wa .^ p.exp_ea;
% capital fraction of egg production
v.capfrac = sum(v.E - v.Einc) ./ sum(v.E);

% check for starvation
isstarving = v.R < -p.rstarv .* v.W;
isalive = isalive & cumsum(isstarving)==0;
isineggprod = isineggprod & isalive;
v.E(~isalive) = 0;


% mortality and survivorship: N(t) ---------------------------------------------
v.m = zeros(NT,NC,NS);
v.m(isalive) = p.m0 .* qg(isalive) .* v.a(isalive) .* v.W(isalive).^(p.theta-1);
	% mort. rate at T, size
v.lnN = cumsum(-v.m) .* dt;
% calculate adult recruitment (= recruitment at the moment egg prod begins)
lnNa = v.lnN;
lnNa(~isineggprod) = nan;
v.Na = exp(max(lnNa));


% contributions to fitness at each t -------------------------------------------
v.dF1 = real(v.E .* exp(v.lnN) .* dt) ./ v.We_theo;
v.dF1(isnan(v.dF1)) = 0;
v.F1 = sum(v.dF1);
v.tEcen = sum(v.t .* v.dF1) ./ v.F1; % center of mass of E*N
	% tEcen - t0 is generation length: similar to, but more accurate than,
	% dtegg - t0
	
	
% two-generation fitness -------------------------------------------------------
F1expected = max(v.F1,[],3); % expected LEP for each t0, assuming that the 
						     % offspring will take the optimal strategy
F1ex_ = interp1(v.t0,F1expected,v.t(:,1)); % expected LEP for offspring produced
										   % on a given timestep
F1ex_(isnan(F1ex_)) = 0;
F1ex_ = repmat(F1ex_(:),[1 NC NS]);
dF2 = v.dF1 .* F1ex_; % contribution to two-generation fitness at each (t,t0,s)
v.F2 = sum(dF2); % two-generation fitness at each (t0,s)
% could iterate like this: F2expected = max(v.F2,[],3); ...
% note that (n-generation fitness)^(1/n) is the per-generation fitness, and
% (n-generation fitness)^(365/n/(tEcen-t0)) - 1 is the per-year growth rate
		    
	
% clean up ---------------------------------------------------------------------
if strcmpi(whatToSave,'scalars only')
	v = keepScalars(v);
elseif strcmpi(whatToSave,'scalars and fitness')
	dF1 = v.dF1;
	v = keepScalars(v);
	v.dF1 = dF1;
else
	v.a(~isalive) = nan;
	v.D(~isalive) = nan;
	v.W(~isalive) = nan;
	v.R(~isalive) = nan;
	v.lnN(~isalive) = -Inf;
	v.E(~isalive) = 0;
end
