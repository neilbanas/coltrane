function v = coltrane_integrate(forcing,p,t0);

% v = coltrane_integrate(forcing,p,t0);
%
% calculates a(t), D(t), W(t), R(t), N(t), E(t), and dF1 for a set of spawning
% dates t0 and a fixed timing strategy (s), along with a variety of summary metrics.
% 
% The three fields of s (tdia_exit, tdia_enter, and dtegg) should be inside p as
% additional scalar parameters.


% copy over forcing time series and put them in the right shape
v = forcing;
v.t0 = t0(:)';
NC = size(v.t0,2); % number of cohorts (t0 values)
NT = size(v.t,1); % number of timesteps
fields = fieldnames(forcing);
for k=1:length(fields)
	v.(fields{k}) = repmat(forcing.(fields{k})(:),[1 NC]);
end

% timebase
v.yday = yearday(v.t); % make sure this is filled in and consistent
dt = v.t(2) - v.t(1);  % determine simulation timestep from the forcing time series

% define output variables that exist for all cases, even if there are no valid
% solutions to the model integration
v.F1 = zeros(1,NC);
v.level = zeros(1,NC);
	% level 0 = logical inconsistency (between development and dtegg)
 	% level 1 = successful diapause
 	% level 2 = reaches adulthood
 	% and past there, fitness provides classifications
 	

% activity: a(t) ------------------------------------------------------
v.a = double(~(v.yday >= p.tdia_enter | v.yday <= p.tdia_exit));
	% assumes that tdia_enter > tdia_exit
isalive = v.t >= repmat(v.t0,[NT 1]); % born yet?
if p.requireActiveSpawning
	% eliminate t0 values fall during diapause
	t0_yday = yearday(v.t0);
	activeSpawning = ~(t0_yday >= p.tdia_enter | t0_yday <= p.tdia_exit);
	isalive = isalive & repmat(activeSpawning,[NT 1]); 
	% likewise t0 values in which t0+dtegg falls during diapause
	% (this eliminates some redundancy)
	ta_yday = yearday(v.t0 + p.dtegg);
	activeNextSpawning = ~(ta_yday >= p.tdia_enter | ta_yday <= p.tdia_exit);
	isalive = isalive & repmat(activeNextSpawning,[NT 1]); 
end

% temperature response factors: qd, qg  --------------------------------
v.temp = v.a .* v.T0 + (1-v.a) .* v.Td;
qd = zeros(NT,NC);
qd(isalive) = (p.Q10d^0.1) .^ v.temp(isalive);
qg = zeros(NT,NC);
qg(isalive) = (p.Q10g^0.1) .^ v.temp(isalive);


% prey saturation ------------------------------------------------------
v = preySaturation(v,p);
if size(v.sat,2)==1
	v.sat = repmat(v.sat(:),[1 NC]);
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
isineggprod = (v.t >= v.t0 + p.dtegg);

% check D(t) for consistency with dtegg
toolate = any(v.D<1 & isineggprod);
v.D(:,toolate) = nan;
if all(toolate)
	% if there are no good solutions at all, don't bother with the rest of the model
	return
end

% calculate D in middle of first winter, just as a diagnostic
[yr0,~] = datevec(v.t0);
first31dec = repmat(reshape(datenum(yr0,12,31), [1 NC]), [NT 1]);
is31dec = abs(v.t-first31dec) == ...
			  repmat(min(abs(v.t-first31dec)),[NT 1]);
is31dec = is31dec & cumsum(is31dec)==1;
v.D_winter = reshape(v.D(is31dec),[1 NC]);

% flag time points at which the animal is in diapause but at a 
% diapause-incapable stage, and mark these cases as dead
isactive = isalive & (v.a==1 | ~isfeeding);
hasbeenactive = (cumsum(isactive) >= 1);
isfailingtodiapause = isalive & hasbeenactive & v.a==0 & v.D < p.Ddia;
hasfailedtodiapause = (cumsum(isfailingtodiapause)>1);
v.D(hasfailedtodiapause) = nan;
v.level(~any(hasfailedtodiapause)) = 1; % level 1 = successful diapause

isalive = isalive & ~isnan(v.D);
isfeeding = isfeeding & isalive;
isineggprod = isineggprod & isalive;

% date on which D=1 is reached (another diagnostic)
tD1 = v.t;
tD1(v.D<1 | ~isfinite(v.D)) = nan;
v.tD1 = min(tD1);
v.level(any(v.D==1)) = 2; % level 2 = reaches adulthood

% pick a guess at adult size based on mean temperature in the forcing,
% and pick a corresponding guess at egg size based on this.
% in v1.0 these were called p.Wa0 and p.We0.
T_nominal = mean(v.temp);
qoverq = (p.Q10g./p.Q10d).^(T_nominal./10);
co = (1-p.theta) .* p.GGE_nominal .* (1-p.Df) .* qoverq;
v.Wa_theo = (co .* p.I0 ./ p.u0) .^ (1./(1-p.theta));
v.We_theo = p.r_ea .* v.Wa_theo .^ p.exp_ea;


% energy gain, growth, egg production: G(t), W(t), R(t), E(t) ------------------
v.G = zeros(NT,NC);
v.W = repmat(v.We_theo,[NT 1]);
v.R = zeros(NT,NC);
v.Einc = zeros(NT,NC);
v.E = zeros(NT,NC);
astar = p.rb + (1-p.rb) .* v.a;

for n=1:NT-1
	f = isfeeding(n,:);
	e = isineggprod(n,:);
	% net gain
	Imax_nf = qg(n,f) .* p.I0 .* v.W(n,f).^(p.theta-1);
	I_nf = v.a(n,f) .* p.r_assim .* v.sat(n,f) .* Imax_nf;
	M_nf = p.rm .* astar(n,f) .* Imax_nf;
	v.G(n,f) = I_nf - M_nf;
    GWdt = v.G(n,:) .* v.W(n,:) .* dt;
	% allocation to growth
	dW = GWdt;
	dW(dW>0 & e) = 0; % definitely
	if ~p.allowGainAfterD1
		dW(dW>0 & v.D(n,:)==1) = 0;
		% positive gain is simply lost between D=1 and the start of egg prod, as if the
		% animals are wishing they were in diapause. This might be too strict, but 
		% otherwise there's no upper limit on the buildup of reserves.
	end
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
last = ~isineggprod(1:end-1,:) & isineggprod(2:end,:);
v.Wa = max(v.W(1:end-1,:) .* last);
v.Ra = max(v.R(1:end-1,:) .* last);
% update the estimate of We
v.We = p.r_ea .* v.Wa .^ p.exp_ea;
% capital fraction of egg production
v.capfrac = sum(v.E - v.Einc) ./ sum(v.E);
% W, R at all stages
stages = {'N6','C1','C2','C3','C4','C5','C6'};
for i=2:length(stages)-1
	Dstart_i = 0.5 .* (stage2D(stages{i-1}) + stage2D(stages{i}));
	Dend_i = 0.5 .* (stage2D(stages{i}) + stage2D(stages{i+1}));
	atstage_i = double(v.D >= Dstart_i & v.D < Dend_i & isalive);
	atstage_i(atstage_i==0) = nan;
	v.(['W_' stages{i}]) = nanmean(v.W .* atstage_i);
	v.(['R_' stages{i}]) = nanmean(v.R .* atstage_i);
end
% W, R for winter late stages
% (this is not a general definition of winter, but calculating it here avoids the
% need to save full time-series output)
D_C5_start = 0.5.*(stage2D('C4') + stage2D('C5'));
iswin = v.yday >= datenum('Nov 1 0000') | v.yday <= datenum('Feb 28 0000');
iswinlatestage = double(v.D >= D_C5_start & isalive & iswin);
iswinlatestage(iswinlatestage==0)=nan;
v.W_C56win = nanmean(v.W.*iswinlatestage);
v.R_C56win = nanmean(v.R.*iswinlatestage);


% check for starvation
isstarving = v.R < -p.rstarv .* v.W;
isalive = isalive & cumsum(isstarving)==0;
isineggprod = isineggprod & isalive;
v.E(~isalive) = 0;


% mortality and survivorship: N(t) ---------------------------------------------
v.m = zeros(NT,NC);
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

% timing metrics ---------------------------------------------------------------
v.tEcen = sum(v.t .* v.dF1) ./ v.F1; % center of mass of E*N
	% tEcen - t0 is generation length: similar to, but more accurate than,
	% dtegg - t0
GWNdt = v.G .* v.W .* exp(v.lnN) .* dt;
GWNdt(isnan(GWNdt)) = 0;
GWNdt = max(0, GWNdt);
v.tGain = sum(v.t .* GWNdt) ./ sum(GWNdt); % center of mass of gain G*W*N
mWNdt = v.m .* v.W .* exp(v.lnN) .* dt;
mWNdt(isnan(mWNdt)) = 0;
v.tYield = sum(v.t .* mWNdt) ./ sum(mWNdt); % center of mass of yield m*W*N
mRNdt = v.m .* v.R .* exp(v.lnN) .* dt;
mRNdt(isnan(mRNdt)) = 0;
v.tYieldR = sum(v.t .* mRNdt) ./ sum(mRNdt); % center of mass of lipid yield m*R*N
mWNdt = v.m .* v.W .* exp(v.lnN) .* dt .* (v.D >= D_C5_start);
mWNdt(isnan(mWNdt)) = 0;
v.tYieldC56 = sum(v.t .* mWNdt) ./ sum(mWNdt); % center of mass of C5-6 yield

	    
% scalar metrics from forcing --------------------------------------------------
fields = {'Ptot','T0','Td','sat','ice','satIA','satWC'};
a0 = v.a;
a0(~isfinite(a0)) = 0;
for k = 1:length(fields)
	if isfield(v,fields{k})
		f0 = v.(fields{k});
		f0(~isfinite(f0)) = 0;
		if size(f0,2)==1, f0 = repmat(f0,[1 NC 1]); end
		v.([fields{k} '_avg']) = nanmean(f0);
		v.([fields{k} '_active']) = sum(f0.*a0) ./ sum(a0);
		isgrowing = isalive & ~isineggprod;
		v.([fields{k} '_growing']) = sum(f0.*isgrowing) ./ sum(isgrowing); 
	end
end	    
	
	
% clean up ---------------------------------------------------------------------
% blank out nonliving portions of the time series
v.a(~isalive) = nan;
v.D(~isalive) = nan;
v.W(~isalive) = nan;
v.R(~isalive) = nan;
v.lnN(~isalive) = -Inf;
v.E(~isalive) = 0;
