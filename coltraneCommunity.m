function coltraneCommunity(outfile,forcing,p0,traits);

% runs Coltrane for a single forcing time series but some range of
% one or more traits, constructing a population for each element of the
% matched fields in _traits_. Each population has a range of spawning dates t0
% and timing strategy tdia_enter, tdia_exit, dtegg.

traitFields = fieldnames(traits);
Ntr = numel(traits.(traitFields{1}));
disp(['varying the traits...']);
disp(traitFields)
disp(['... in ' num2str(Ntr) ' combinations.']);

for i=1:Ntr
	disp(['trait combination ' num2str(i) ' / ' num2str(Ntr)]);

	% parameter set for trait combination i
	p = p0;
	for k=1:length(traitFields)
		p.(traitFields{k}) = traits.(traitFields{k})(i);
	end
	
	% run one population
	vi = coltranePopulation(forcing,p,'scalars and fitness');
	
	% we save the full fitness landscape in coltrane_integrate.m so that 
	% coltranePopulation.m can calculate the 2-generation fitness--
	% but after that it's unnecessary:
	vi = rmfield(vi,{'dF1','t'});
	
	% A.* = all cohorts/strategies
%	fields = fieldnames(vi); % to save everything
	fields = {'F1','F2','Wa','capfrac','tEcen','D_winter','t0','level'};
	for k=1:length(fields)
		vifk = squeeze(vi.(fields{k}));
		vifk = reshape(vifk,[1 size(vifk)]);
		if i==1
			A.(fields{k}) = repmat(nan.*vifk,[Ntr 1 1]);
		end
		A.(fields{k})(i,:,:) = vifk;
	end
end



% clean up A and save -----------------
try
fields = fieldnames(A);
	for k=1:length(fields)
		A.(fields{k}) = reshape(A.(fields{k}),[size(traits.(traitFields{1})) ... 	
											size(A.(fields{k}),2) size(A.(fields{k}),3)]);
	end
catch
	warning('reshaping problem with A.');
end
save('-v7.3',allFile, '-struct', 'A');
