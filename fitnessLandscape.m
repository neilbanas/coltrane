function fl = coltrane_fitnessLandscape(dF1,t,t0,s);

fl.dF1 = dF1;
fl.t = t;

[NT,NC,NS] = size(dF1);
	% NT should match length of t
	% NC should match length of t0
	% NS should match length of fields of s

fields = fieldnames(s);
sz = size(s.(fields{1})); % prod(sz) = NS
fl.t0 = repmat(t0(:),[1 sz]);
for k=1:length(fields)
	fl.(fields{k}) = repmat(reshape(s.(fields{k}),[1 sz]),[NC 1]);
end

fl.F1 = sum(dF1); % lifetime egg production for each (t0,s)
fl.F1expected = max(fl.F1,[],3); % expected LEP for each t0, assuming that the 
						   % offspring will take the optimal strategy
F1ex_ = interp1(t0,fl.F1expected,t);
F1ex_(isnan(F1ex_)) = 0;
F1ex_ = repmat(F1ex_(:),[1 NC NS]);
fl.dF2 = dF1 .* F1ex_; % contribution to two-generation fitness at each (t,t0,s)
fl.F2 = sum(fl.dF2); % two-generation fitness at each (t0,s)
fl.F2expected = max(fl.F2,[],3); % at each t0, assuming that the offspring's
						   % offspring take the optimal strategy
% could iterate like this...

