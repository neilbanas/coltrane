function coltraneCommunity(outfile,forcing,p0,traits,doFiltering);

% coltraneCommunity(outfile,forcing,p0,traits,doFiltering);
%
% runs Coltrane for a single forcing time series but some range of
% one or more traits, constructing a population for each element of the
% matched fields in _traits_. Each population has a range of spawning dates t0
% and timing strategy s = (tdia_enter, tdia_exit, dtegg).
%
% if doFiltering = 'no filtering' or 0, returns every combination of traits x t0 x s
% otherwise, ***not written yet, wasn't there something about sequences in dtegg
% *** also, only return t0 in the first year, right?

traitFields = fieldnames(traits);
Ntr = numel(traits.(traitFields{1}));
disp(['varying the traits...']);
disp(traitFields)
disp(['... in ' num2str(Ntr) ' combinations.']);


% run each trait combination ------------------

for i=1:Ntr
	disp(['trait combination ' num2str(i) ' / ' num2str(Ntr)]);

	% parameter set for trait combination i
	p = p0;
	for k=1:length(traitFields)
		p.(traitFields{k}) = traits.(traitFields{k})(i);
	end
	
	% run one population
	pop_i = coltranePopulation(forcing,p);
	
	% add results to one big structure
	fields = fieldnames(pop_i); % to save everything
%	fields = {'F1','F2','Wa','capfrac','tEcen','D_winter','t0','level'}; % highlights only
	for k=1:length(fields)
		if ~isfield(pop_i,fields{k}), continue; end
		pop_i_var_k = squeeze(pop_i.(fields{k}));
		pop_i_var_k = reshape(pop_i_var_k,[1 size(pop_i_var_k)]);
		if i==1 % initialise one field of the output structure _comm_
			comm.(fields{k}) = repmat(nan.*pop_i_var_k,[Ntr 1 1]);
		end
		comm.(fields{k})(i,:,:) = pop_i_var_k; % save one slice of _comm_
	end
	
end


% clean up output and save -----------------

% add traits to main output structure
fields = fieldnames(traits);
sz = size(comm.F1);
	for k=1:length(fields)
	comm.(fields{k}) = repmat(traits.(fields{k})(:),[1 sz(2:end)]);

end

save('-v7.3',outfile, 'comm','forcing','p0','traits');
