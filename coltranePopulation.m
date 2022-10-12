function v = coltranePopulation(forcing,p,whatToSave);

% v = coltranePopulation(forcing, p, ...
%           'scalars only' | 'scalars and fitness' | 'everything');
%
% v2.0 of the Coltrane model. This has diverged significantly from the
% Coltrane 1.0 model in Banas et al., Front. Mar. Res., 2016.
%
% forcing is a structure specifying a single time series of forcing and 
% 		ancillary variables.
% p is a structure containing the internal model parameters.
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


if nargin < 3, whatToSave = 'scalars only'; end

fields = fieldnames(forcing);
for i=1:length(fields)
	forcing.(fields{i}) = forcing.(fields{i})(:);
end
NT = size(forcing.t,1); % # timesteps
[t0,s] = timingCombinations(forcing,p);
strategyFields = fieldnames(s);
NC = length(t0);
NS = prod(size(s.(strategyFields{1})));


% run one strategy at a time
out = cell(1,NS);
disp([num2str(NT) ' timesteps x ' num2str(NC) ' cohorts x ' num2str(NS) ' strategies']);
for i = 1:NS
	pii = addStrategyToParams(p,s,i);
	out{i} = coltrane_integrate(forcing,pii,t0,whatToSave);
end


% rearrange into something more useful. (_parfor_ prefers that we don't do all this
% fancy indexing inside the main loop over strategies 1:NS.)
for i = 1:NS
	v.level(:,i) = out{i}.level(:);
end
f = find(any(v.level>0)); % strategies with any complete integrations
if isempty(f)
	% if there aren't any, return only level (and F1, F2 = 0)
	return;
end
fields = fieldnames(out{f(1)});
for k = 1:length(fields)
	v.(fields{k}) = repmat(nan,[size(out{f(1)}.(fields{k})) NS]);
	for i=1:NS
		if isfield(out{i},fields{k}) && ~isempty(out{i}.(fields{k}))
			v.(fields{k})(:,:,i) = out{i}.(fields{k});
		end
	end
end


% two-generation fitness
if isfield(v,'dF1')
	F1expected = max(v.F1,[],3); % expected LEP for each t0, assuming that the 
								 % offspring will take the optimal strategy
	F1ex_ = interp1(v.t0(:,:,1),F1expected,v.t(:,1));
		% expected LEP for offspring produced on a given timestep
	F1ex_(isnan(F1ex_)) = 0;
	F1ex_ = repmat(F1ex_(:),[1 NC NS]);
	dF2 = v.dF1 .* F1ex_; % contribution to two-generation fitness at each (t,t0,s)
	v.F2 = sum(dF2); % two-generation fitness at each (t0,s)
	% could iterate like this: F2expected = max(v.F2,[],3); ...
	% 
	% note that (n-generation fitness)^(1/n) is the per-generation fitness, and
	% (n-generation fitness)^(365/n/(tEcen-t0)) - 1 is the per-year growth rate
	v.F1 = squeeze(v.F1);
	v.F2 = squeeze(v.F2);
end


for k = 1:length(fields)
	v.(fields{k}) = squeeze(v.(fields{k}));
end

