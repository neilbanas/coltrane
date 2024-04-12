function comm_out = filterCommunity(comm);

% filter output from coltraneCommunity down to the small fraction of cases
% that are useful. This would preclude some types of analysis, like examining
% "level" to see why failed cases failed.

% keep only first year of t0 values
isFirstYear = comm.t0 < comm.t0(1) + 365;
% keep only F1 > 0
nonzeroF1 = comm.F1>0;
% two subsets of what remains:
% - best F1 for each trait x t0
bestFor_t0 = comm.F1==max(comm.F1,[],3);
% - other cases that are within 0.1x of the global max F1.
%	this is a quick replacement for identifying local maxima along the
%   dtegg axis
okRelativeToGlobalMax = comm.F1>0.1*max(comm.F1(:));

f = find(isFirstYear & nonzeroF1 & (bestFor_t0 | okRelativeToGlobalMax));

fields = fieldnames(comm);
for k = 1:length(fields)
	comm_out.(fields{k}) = comm.(fields{k})(f);
end
