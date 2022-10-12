function pii = addStrategyToParams(p,s,ind)

% pii = addStrategyToParams(p,s,ind);

pii = p;
fields = fieldnames(s);
for k=1:length(fields)
	pii.(fields{k}) = s.(fields{k})(ind);
end