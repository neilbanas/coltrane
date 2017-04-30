function [cohort,pop,routine] = coltranePostproc_ER(v,p,forcing);

% postprocessing for the Coltrane egg/reserve model
% (see Banas et al., Front. Mar. Sci., 2016)


% metrics for each cohort ------------------------------------------------------
% (cohort.* are raw metrics for each spawning date; pop.* are metrics for the
% population, once a stable annual cycle of spawning has been determined)

% actual adult size
Na = exp(v.lnN);
SR = v.S+v.R;
Na(v.D<1 | isnan(SR)) = 0;
SR(isnan(SR)) = 0;
cohort.Wa = sum(SR.*Na)./sum(Na); % abundance-weighted average of adult states
cohort.We = p.r_ea .* cohort.Wa.^p.exp_ea;
	% updated egg size: if Wa = v.Wa_theo, We = v.We_theo

% lifetime egg production: expected eggs per egg
cohort.LEP_by_egg = sum((v.Einc+v.Ecap).*exp(v.lnN))./cohort.We.*p.dt;
cohort.LEP_by_egg(isnan(cohort.LEP_by_egg)) = 0;
% note: a simple way to take into account internal life-history mismatch without
% the eigenvalue calculation below would be a two-generation calculation:
% instead of LEP = integral of E * N, LEP2 = integral of E * N * LEP

% terms in the lifetime energy budget
% yield to predators = total gain - egg prod beyond replacement
de = (v.Fin-v.Fmet).*v.S.*exp(v.lnN);
de(isnan(de)) = 0;
cohort.ener_gain = p.dt .* sum(de);
cohort.ener_egg = (cohort.LEP_by_egg - 1).*cohort.We;
de = v.m.*(v.S+v.R).*exp(v.lnN);
de(isnan(de)) = 0;
cohort.ener_yield = p.dt .* sum(de); % for checking
de = v.m.*v.R.*exp(v.lnN);
de(isnan(de)) = 0;
cohort.ener_yield_R =  p.dt .* sum(de);
cohort.rhoY = cohort.ener_yield_R ./ cohort.ener_yield;
	% yield (trophic transfer) of R as fraction of total yield
cohort.ener_starv = cohort.ener_gain - cohort.ener_egg - cohort.ener_yield;
	% the residual of the three terms is (R+S)N at the moment of starvation--
	% accumulation of energy gain in the cohort itself

% other metrics by cohort
cohort.develtime = sum(v.D < 1) .* p.dt;
D1 = v.D;
D1(v.a~=0) = nan;
cohort.Ddiamin = min(D1);
rho = max(0,v.R) ./ (v.R + v.S);
WN = (v.R+v.S).*exp(v.lnN);
WN(isnan(WN))=0;
rhoWN = rho .* WN;
rhoWN(isnan(rhoWN))=0;
cohort.rhomean = sum(rhoWN) ./ sum(WN);
Pgrow = forcing.P;
Pgrow(v.a < 1 | v.D >=1) = nan;
cohort.satmean = nanmean(Pgrow ./ (p.Ks + Pgrow));
cohort.capfrac = nanmean(v.Ecap.*exp(v.lnN)) ./ ...
				 nanmean((v.Ecap+v.Einc).*exp(v.lnN));
NA = exp(v.lnN);
NA(v.D<1) = 0;
cohort.recruitmentToAdult = max(NA);


% find a stable cycle of egg production ----------------------------------------

% V0 ~ EN/We dt; sum over first dimension (time) is LEP
N = length(v.t);
NC = length(v.tspawn);
pop.V0 = (v.Einc+v.Ecap).*exp(v.lnN)./repmat(cohort.We,[N 1]).*p.dt;
pop.V0(isnan(pop.V0)) = 0;
% chunk it into an annual cycle (by tspawn, ignoring which year a given bit of 
% spawning happens). Make sure p.dt is a divisor of p.dt_tspawn, and both fit
% evenly into a 365-day year!
pop.V = zeros(NC,NC);
i0 = repmat((1:(p.dt_tspawn/p.dt))',[1 p.Nyears]);
i0 = i0 + repmat((NC-1).*size(i0,1).*(0:p.Nyears-1),[size(i0,1) 1]);
i0 = i0(:);
for j=1:NC
	pop.V(j,:) = pop.V(j,:) + sum(pop.V0(i0 + (j-1)*(p.dt_tspawn/p.dt), :));
end
cohort.genlength = sum((forcing.t - repmat(v.tspawn(:)',[N 1])) .* pop.V0) ...
				./ sum(pop.V0); % time from tspawn to mean egg prod date

% pop.V is a transition matrix for generating an annual time series of eggs 
% spawned from the previous year of eggs spawned. The first eigenvector of pop.V 
% is thus a stable annual cycle of egg production.
% n(t)_next = V' * n(t)'
[evecs,evals] = eig(pop.V');
eval1 = max(real(diag(evals)));
if eval1 > 0
	pop.lam = eval1 - 1; % first eigenvalue is 1 + population growth rate lambda
	f = find(real(diag(evals))==eval1);
	pop.n = evecs(:,f(1))';
	pop.n(isnan(pop.n)) = 0;
	pop.n = pop.n ./ sum(pop.n);
		% normalize so that pop.n is _relative_ distribution of spawning dates--
		% a set of weights, and therefore proportional to p.dt_tspawn
else
	pop.lam = nan;
	pop.n = zeros(size(v.tspawn));
end
	
	
% population-level metrics -----------------------------------------------------
fields = fieldnames(cohort);
for i=1:length(fields)
	F = cohort.(fields{i});
	F(isnan(F))=0;
	pop.(fields{i}) = sum(F.*pop.n) ./ sum(pop.n);
		% previous versions weighted by n * LEP instead of n, which is perhaps
		% overthinking it
end
tt = [v.tspawn(cumsum(pop.n)>0.1.*sum(pop.n)) nan];
pop.tE10 = tt(1); % yearday when 1st 10% of egg prod completed
tt = [v.tspawn(cumsum(pop.n)>0.5.*sum(pop.n)) nan];
pop.tE50 = tt(1);
tt = [v.tspawn(cumsum(pop.n)>0.9.*sum(pop.n)) nan];
pop.tE90 = tt(1);
pop.LEP_by_adult = pop.LEP_by_egg ./ pop.recruitmentToAdult;
	% lifetime egg production per successfully recruited adult
	% (per adult female might be twice this)
pop.Wa_theo = v.Wa_theo;
pop.We_theo = v.We_theo;
	 % preserve the initial guesses at adult and egg size (see beginning of
	 % coltraneModel.m, before the main loop)
	 

% annual routine ---------------------------------------------------------------
% abundance-weights for adding up time series to make a population.
% Nn = relative spawning intensity n * survivorship of each case N, taking into 
% account carryover from previous years (this last point is what makes the code 
% a bit complicated)
kk = 1:(p.dt_tspawn/p.dt):(365/p.dt);
rsh = [365/p.dt p.Nyears+1];
Nn0 = exp(v.lnN) .* repmat(pop.n,[N 1]);
Nntot = sum(sum(reshape(Nn0(1:end-1,:)',[NC rsh]),3),1);
Nn = [Nn0(1:end-1,:) ./ repmat(Nntot(:),[p.Nyears+1 NC]); zeros(1,NC)];
% another version, for diapause-capable stages only
Nndia0 = exp(v.lnN) .* repmat(pop.n,[N 1]);
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
routine.stageFracN = stageFracN0(kk,:);
routine.stageFracW = stageFracW0(kk,:);
% other population time series
afrac0 = sum((v.a==1) .* Nndia, 2);
afrac0 = sum(reshape(afrac0(1:end-1),rsh),2);
routine.afrac = afrac0(kk);
	% fraction of diapause-capable individuals that are in diapause
	
pop.afrac = mean(routine.afrac);
	% overall fraction of days * diapause-capable individuals that are active

routine.t = v.tspawn(1:end-1);
routine.relativeEggProd = pop.n(1:end-1);
	% relative egg production is just n, but it's copied over here for clarity