function v = filterCohorts(v_in,criteriaSet);

% v_out = filterCohorts(v_in,criteriaSet);
%
% applies some set of critera to coltrane output and returns a smaller selection
% of cases. If the fields of v_in are size [NT NP NC], then fields of v_out are
% size [NT NF], where NF <= NP*NC.

v = v_in;
[NT,NP,NC] = size(v.t);


if strcmpi(criteriaSet,'bcc')
	% some filters appropriate for C. glacialis/marshallae in the
	% Bering/Chukchi system
	startOfC5 = (stage2D('C4')+stage2D('C5'))/2;
	develTests = v.D_winter > startOfC5 ... % overwinters as C5
			   & v.numWinters <= 1 ...		% overwinters only once
			   & v.activeSpawning;			% doesn't spawn during diapause
	energeticsTests = v.Wa < 252 + 2*65 ... % mean +/- 2 s.d.,
					& v.Wa > 252 - 2*65 ... % 		Campbell et al. 2016
					& v.starv_stress < 0.5; % just a guess
	popDynTests = v.F > 1 ...      % replacement rate
				& v.capfrac < 0.5; % predominantly income breeding
	ind = find(develTests & energeticsTests & popDynTests);
else
	warning(['don''t recognise ''' criteriaSet ''', so ignoring it.']);
	ind = 1:(NP*NC);
end


% do the filtering
fields = fieldnames(v);
for i=1:length(fields)
	v.(fields{i}) = v.(fields{i})(:,ind);
end

