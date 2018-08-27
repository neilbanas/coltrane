function v = coltrane_oneMaturationDate(v0,p,tadult,forceCompleteSolution);

% v = coltrane_oneMaturationDate(v0,p,tadult);
%
% calculates a(t), D(t), W(t), N(t), E(t), and LEP1 for a single value of
% tadult.
%
% see coltraneModel.m for fuller notes

if nargin < 4, forceCompleteSolution = 0; end

v = v0;
% assumes forcing time series have been copied over and t0, tdia_enter,
% tdia_exit have been filled in
dt = v.t(2) - v.t(1);
[NT,NC,NDx,NDn] = size(v.t);

% activity: a(t) ------------------------------------------------------
v.a = double(~(v.yday >= repmat(v.tdia_enter,[NT 1 1 1]) | ...
			   v.yday <= repmat(v.tdia_exit,[NT 1 1 1])));
	% assumes that tdia_enter > tdia_exit
								
% temperature response factors: qd, qg  --------------------------------
v.temp = v.a .* v.T0 + (1-v.a) .* v.Td;
qd = (p.Q10d^0.1) .^ v.temp;
qg = (p.Q10g^0.1) .^ v.temp;

% prey saturation ------------------------------------------------------
v = preySaturation(v,p);

% development: D(t) ----------------------------------------------------
isalive = v.t >= repmat(v0.t0,[NT 1 1 1]); % spawned yet?
dDdt = isalive .* p.u0 .* qd; % nonfeeding formula
v.D = cumsum(dDdt) .* dt;
isfeeding = v.D >= p.Df;
dDdt_feeding = isalive .* p.u0 .* qd .* v.a .* v.sat;
dDdt(isfeeding) = dDdt_feeding(isfeeding); % full formula 
v.D = cumsum(dDdt) .* dt;
v.D(v.D>1) = 1;

% check D(t) against tadult
isadult = v.t >= tadult;
toolate = any(v.D<1 & isadult);
if all(toolate) & ~forceCompleteSolution
	v.D = [];
	return;
	% if there are no good solutions at all, don't bother with the rest of the
	% model. isempty(v.D) is the test of whether the model completed
else
	v.D(:,any(v.D<1 & isadult)) = nan;
	% the easy way to harmonise them is to blank out the cases where
	% tadult is too early to be a possible date of maturation given D(t) 
end

% pick a guess at adult size based on mean temperature in the forcing,
% and pick a corresponding guess at egg size based on this.
% in v1.0 these were called p.Wa0 and p.We0.
T_nominal = mean(v.temp);
qoverq = (p.Q10g./p.Q10d).^(T_nominal./10);
co = (1-p.theta) .* p.GGE_nominal .* (1-p.Df) .* qoverq;
v.Wa_theo = (co .* p.I0 ./ p.u0) .^ (1./(1-p.theta));
v.We_theo = p.r_ea .* v.Wa_theo .^ p.exp_ea;

% juvenile growth: W(t) ------------------------------------------------
% approximate curve of W(D), valid for a=1
% this is used for the allometry in the calculation of net gain
isgrowing = v.D >= p.Df & ~isadult;
satmean = sum(v.sat .* isgrowing) ./ sum(isgrowing);
c = (1-p.theta) .* qg ./ qd .* p.I0 ./ p.u0 .* ...
	(p.r_assim - p.rm ./ repmat(satmean,[NT 1 1 1]));
c = max(c,0);
Wapprox = zeros(size(v.D));
We_theo_4d = repmat(v.We_theo,[NT 1 1 1]);
Wapprox(isgrowing) = ((We_theo_4d(isgrowing).^(1-p.theta)) + ...
	c(isgrowing) .* (v.D(isgrowing) - p.Df)) .^ (1/(1-p.theta));
% net gain over development (here stored as gain * W)
ImaxW = zeros(size(Wapprox));
ImaxW(isgrowing) = ...
	qg(isgrowing) .* p.I0 .* Wapprox(isgrowing) .^ p.theta;
astar = p.rb + (1-p.rb) .* v.a;
GW = zeros(size(v.D));
GW(isgrowing) = ImaxW(isgrowing) .* ...
	(v.a(isgrowing) .* p.r_assim .* v.sat(isgrowing) - ...
	p.rm .* astar(isgrowing));
% integrate to get the correct W over development
v.W = We_theo_4d + cumsum(GW) .* dt;
v.W = max(0,v.W);
% adult size Wa
last = v.D(1:end-1,:,:,:) < 1 & v.D(2:end,:,:,:)==1;
v.Wa = max(v.W(1:end-1,:,:,:) .* last);
Wa_4d = repmat(v.Wa,[NT 1 1 1]);
v.W(isadult) = Wa_4d(isadult);
% calculate net gain during adulthood (used later to calculate egg 
% production) and fill in juvenile ingestion, metabolism, and net gain
% while we're at it
Imax = qg .* p.I0 .* v.W.^(p.theta-1);
Imax(~isfeeding) = 0;
Imax(v.W < We_theo_4d) = 0;
v.I = v.a .* p.r_assim .* v.sat .* Imax;
	% this was r_assim * I, as opposed to I, in Coltrane 1.0
v.M = p.rm .* astar .* Imax;
v.G = v.I - v.M;
v.Wrel = v.W ./ cummax(v.W);
	% size relative to largest size previously attained

% mortality and survivorship: N(t) -------------------------------------
v.m = p.m0 .* qg .* v.a .* v.W.^(p.theta-1); % mort. rate at T, size
v.m(~isalive) = 0;
v.lnN = cumsum(-v.m) .* dt;
% calculate adult recruitment
lnNa = v.lnN;
lnNa(~isadult) = nan;
v.Na = exp(max(lnNa));

% egg production -------------------------------------------------------
% income egg prod = GW during the adult period
v.Einc = max(0, v.G .* v.W);
v.Einc(~isadult) = 0;
% total capital egg prod is a sensible fraction of body size at the end of the
% juvenile period, and then we spread this out over an interval dt_spawn
% starting at t = tadult
Ecap = zeros(size(v.Einc));
sum_Ecap_dt = p.capitalEfficiency .* v.Wa;	
ncap = v.t(:,1) >= tadult & v.t(:,1) < tadult + p.dt_spawn;
num_ncap = sum(ncap);
Ecap(ncap,:,:,:) = repmat(sum_Ecap_dt./num_ncap./dt,[num_ncap 1 1 1]);
v.E = v.Einc + Ecap;
v.capfrac = sum(Ecap) ./ sum(v.E);

% fitness: one-generation calculation ----------------------------------
% ignores timing and internal mismatch
v.LEP1 = sum(v.E .* exp(v.lnN) .* dt) ./ v.We_theo;

% metrics of viability and other classifications -----------------------
[yr0,~] = datevec(v.t0);
yr0 = reshape(yr0,[1 NC NDx NDn]);
[yra,~] = datevec(tadult);
v.numWinters = (yra - yr0);
	% how many winters to reach adulthood (should be <= 1)
first31dec = repmat(reshape(datenum(yr0,12,31), [1 NC NDx NDn]), [NT 1 1 1 1]);
is31dec = abs(v.t-first31dec) == ...
			  repmat(min(abs(v.t-first31dec)),[NT 1 1 1 1]);
v.D_winter = reshape(v.D(is31dec),[1 NC NDx NDn]);
	% D in middle of first winter (should be > params.Ddia)
v.starv_stress = max((1 - v.Wrel)./v.D);
	% starvation stress. % this should be less than some threshhold c
	% which is equivalent to the criterion Wrel >= 1 - c*D
	% (e.g., if c=0.5, then animals are allowed to metabolise half their 
	% body mass at adulthood)
		
% define viability (an additional filter on LEP1)
% v.viable = v.activeSpawning & v.starv_stress < 1;
v.viable = v.starv_stress < 1;
v.LEP1 = v.LEP1 .* v.viable;
v.LEP1(isnan(v.LEP1)) = 0;

% weightings alone the diapause-strategy dimensions --------------------
% relative contributions of each diapause strategy to the total fitness
% of each t0 cohort (at a given tadult)
v.diaStrategyFrac = v.LEP1 ./ ...
	repmat(sum(sum(v.LEP1,3),4), [1 1 NDx NDn]);
v.diaStrategyFrac(isnan(v.diaStrategyFrac)) = 0;
% alternate weightings that just average together the ones that are tied for
% max LEP
maxLEP1 = repmat(max(max(v.LEP1,[],3),[],4),[1 1 NDx NDn]);
nummax = repmat(sum(sum(maxLEP1,3),4),[1 1 NDx NDn]);
v.diaStrategyFrac_opt = double(v.LEP1 == maxLEP1) ./ nummax;
v.diaStrategyFrac_opt(isnan(v.diaStrategyFrac_opt)) = 0;

% clean up -------------------------------------------------------------
fields = fieldnames(v);
for i=1:length(fields)
	if size(v.(fields{i}),1) == NT && ~strcmpi(fields{i},'t')
		v.(fields{i})(~isalive) = nan;
	end
end

