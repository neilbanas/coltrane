function [pop,popts] = coltranePopulation(forcing,p);

% pop = coltranePopulation(forcing, p);
% [pop,popts] = coltranePopulation(forcing, p);
%
% v2.0 of the Coltrane model. This has diverged significantly from the
% Coltrane 1.0 model in Banas et al., Front. Mar. Res., 2016.
%
% forcing is a structure specifying a single time series of forcing and 
% 	ancillary variables.
% p is a structure containing the internal model parameters.
%
% pop contains scalar summaries of what happened in each cohort/strategy combo
% popts contains full time series of state variables, dF1, and so on.
%   To save space, forcing is omitted from these.
%
% relationship with the published Coltrane model (Front. Mar. Res. 2016):
% * the phi model is hereby abandoned
% * R,S have been replaced by R and W = S + R. Allometric formulas using S
%   in the paper now use W.
% * the state variables have been largely separated so that not all parts
% 		require iteration in time. This makes the model cleanly 
% 		hierarchical, so that one can derive predictions about
%				1) development alone (D),
% 				2) then size and time evolution in surviving cohorts (D,W,R),
% 				3) then mortality, survivorship, and population dynamics (D,W,R,N,E).
%		This will also make it possible for N to be density-dependent in a 
%		future version.
% * The myopic criterion for diapause has been replaced by a matrix of entry
%		and exit dates, which are analyzed in a brute-force way parallel
% 		to spawning date.
% * tegg has been replaced by dtegg, which is similar to (tegg - t0).
%
%	t			t0			tdia_exit	tdia_enter	dtegg
%	timestep	spawn date	exit date	entry date	egg prod date
%	(calendar)	(calendar)	(yearday)	(yearday)	(relative to t0)
%
% the last three of these are folded into a single strategy vector s.
%
% *** in process Oct 2022: moving to a cleaner hierarchy,
% integrate -> Population -> Community 
% *** to do: save time series in their own variable


retainTimeSeries = (nargout >= 2);
% if popts is requested as an output, we save full time series (which get really big),
% otherwise discard them and save only summaries

NT = size(forcing.t(:),1); % # timesteps
[t0,s] = timingCombinations(forcing,p);
strategyFields = fieldnames(s);
NC = length(t0);
NS = prod(size(s.(strategyFields{1})));


% run one strategy at a time -----------------------------------------------------------
out = cell(1,NS);
disp([num2str(NT) ' timesteps x ' num2str(NC) ' cohorts x ' num2str(NS) ' strategies']);
parfor i = 1:NS
	pii = addStrategyToParams(p,s,i);
	out{i} = coltrane_integrate(forcing,pii,t0);
end


% clean up output ----------------------------------------------------------------------

% figure out which strategies produced successful cases
for i = 1:NS
	v.level(:,i) = out{i}.level(:);
end
f = find(any(v.level>0)); % strategies with any complete integrations
if isempty(f)
	% if there aren't any, return only level (and F1, F2 = 0)
	v.F1 = zeros(size(level));
	v.F2 = zeros(size(level));
	return;
end

% rearrange the output from all the individual runs into a single structure
% (plus a second one for time series variables if we're keeping them)
forcingFields = cat(1,fieldnames(forcing),'yday'); % time series to always omit
fields = fieldnames(out{f(1)});
for k = 1:length(fields)
	example = out{f(1)}.(fields{k}); % our example of this var
	if size(example,1)==1 % if it's summary field, size 1 x something
		pop.(fields{k}) = repmat(nan,[size(example) NS]);
		for i=1:NS
			if isfield(out{i}, fields{k}) && ~isempty(out{i}.(fields{k}))
				pop.(fields{k})(:,:,i) = out{i}.(fields{k});
			end
		end
	elseif strcmpi(fields{k},'dF1') || strcmpi(fields{k},'t') || ...
		(retainTimeSeries && ~ismember(fields{k},forcingFields))
		% if it's a time series, and we're retaining time series, and it isn't in the
		% always-omit list -- or if it's dF1 or t, which we always save since they're 
		% necessary to compute the two-generation fitness:
		popts.(fields{k}) = repmat(nan,[size(example) NS]);
		for i=1:NS
			if isfield(out{i}, fields{k}) && ~isempty(out{i}.(fields{k}))
				popts.(fields{k})(:,:,i) = out{i}.(fields{k});
			end
		end
	end
end


% two-generation fitness ----------------------------------------------------------------
if isfield(popts,'dF1')
	F1expected = max(pop.F1,[],3); % expected LEP for each t0, assuming that the 
								 % offspring will take the optimal strategy
	F1ex_ = interp1(pop.t0(:,:,1),F1expected,popts.t(:,1));
		% expected LEP for offspring produced on a given timestep
	F1ex_(isnan(F1ex_)) = 0;
	F1ex_ = repmat(F1ex_(:),[1 NC NS]);
	dF2 = popts.dF1 .* F1ex_; % contribution to two-generation fitness at each (t,t0,s)
	pop.F2 = sum(dF2); % two-generation fitness at each (t0,s)
	% could iterate like this: F2expected = max(v.F2,[],3); ...
end
% note that (n-generation fitness)^(1/n) is the per-generation fitness, and
% (n-generation fitness)^(365/n/(tEcen-t0)) - 1 is the per-year growth rate



% final cleanup -------------------------------------------------------------------------
% All the fields of pop are size [1 NC NS], which is consistent with the output from
% coltrane_integrate and with the time series in popts (which are [NT NC NS]), but
% annoying to work with. So make them [NC NS].
fields = fieldnames(pop);
for k = 1:length(fields)
	pop.(fields{k}) = squeeze(pop.(fields{k}));
end
% include the strategy fields in pop. (Since coltrane_integrate sees them as model
% parameters, they aren't saved in the output structures manipulated above).
for k = 1:length(strategyFields)
	pop.(strategyFields{k}) = repmat(s.(strategyFields{k})(:)',[NC 1]);
end



% --------------------------------------------------------------------------------------


% copy a single value of the strategy vector s into p as parameters.
% complicated indexing, so parfor prefers that we do it inside a
% standalone function.
function pii = addStrategyToParams(p,s,ind)
pii = p;
fields = fieldnames(s);
for k=1:length(fields)
	pii.(fields{k}) = s.(fields{k})(ind);
end
