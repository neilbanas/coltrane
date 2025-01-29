function comm_out = filterCommunity(comm,severity);

% filter output from coltraneCommunity down to the small fraction of cases
% that are useful. This would preclude some types of analysis, like examining
% "level" to see why failed cases failed.

if nargin < 2, severity = 'light'; end

% keep only first year of t0 values
isFirstYear = comm.t0 < comm.t0(1) + 365;
% keep only F1 > 0
nonzeroF1 = comm.F1>0;

if strcmpi(severity,'heavy')
	% a few subsets of what remains:
	% - best F1 for each trait x t0
	bestF1_by_t0 = comm.F1==max(comm.F1,[],3);
	bestF2_by_t0 = comm.F2==max(comm.F2,[],3);
	% - other cases that are within 0.1x of the global max F1.
	%	this is a quick replacement for identifying local maxima along the
	%   dtegg axis
	okRelativeToGlobalMax = comm.F1>0.1*max(comm.F1(:));
	% better (but unwritten):
	% forget that last thing; instead, find best F1 and F2 by (t0,gl) combination

	f = find(isFirstYear & nonzeroF1 & ...
	     (bestF1_by_t0 | bestF2_by_t0 | okRelativeToGlobalMax));

else % light filtering
	f = find(isFirstYear & nonzeroF1);
end

% ----

% do the filtering
fields = fieldnames(comm);
for k = 1:length(fields)
	comm_out.(fields{k}) = comm.(fields{k})(f);
end