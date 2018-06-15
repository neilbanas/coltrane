function v1 = keepScalars(v);

% v1 = keepScalars(v);
%
% copies over only the fields of v that are length 1 in the first dimension,
% squeezing away that dimension. This is useful for keeping only the diagnostics
% from a coltrane run and discarding the full time series.

fields = fieldnames(v);
for i=1:length(fields)
	if size(v.(fields{i}),1)==1
		v1.(fields{i}) = squeeze(v.(fields{i}));
	end
end
