function si = selectRows(s,ind)

% si = selectRows(s,ind);

fields = fieldnames(s);
for k=1:length(fields)
	si.(fields{k}) = s.(fields{k})(ind);
end