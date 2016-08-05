filename = '/Users/neil/datasets/coltrane_runs/disko/';
Nyears = 5;
paramCases = {'Nyears',Nyears,...
			  'tegg',0:5:365*4,...
			  'u0',0.004:0.0005:0.01,...
			  'Ks',1,...
			  'm0',0.06};
forcingCases = {'dT',0};


% run the runs
if 1
	[cases,stats,p1,forcing1] = ...
		coltraneEnsemble('disko',forcingCases,paramCases,filename);
	save runs_disko_stats stats cases p1 forcing1
else
	load runs_disko_stats stats cases p1 forcing1
end

% identify species analogs
clear spp
spp.fin.W = 140;
spp.gla.W = 370;
spp.hyp.W = 1600;
spp.names = fieldnames(spp);
tol = 0.3;
f = stats.LEPn > 1;
for i=1:length(spp.names)
	sp = spp.names{i};
	spp.(sp).f = f & ...
			   stats.Wa>=spp.(sp).W*(1-tol) & ...
			   stats.Wa<=spp.(sp).W*(1+tol);
	spp.(sp).f1 = find(spp.(sp).f & ...
					   stats.LEPn==max(stats.LEPn(spp.(sp).f)));
	spp.(sp).f1 = spp.(sp).f1(1);

	fname = num2filename(filename,spp.(sp).f1);
	load(fname,'v');
	spp.(sp).stageFracW1 = v.stageFracW;
	spp.(sp).n1 = v.n;
	spp.(sp).afrac1 = v.afrac;
	spp.(sp).v1 = v;
	
	clear stageFracW n afrac
	ff = find(spp.(sp).f);
	for j=1:length(ff)
		fname = num2filename(filename,ff(j));
		load(fname,'v');
		spp.(sp).stageFracW(:,:,j) = v.stageFracW;
		spp.(sp).n(:,j) = v.n(:);	
		spp.(sp).afrac(:,j) = v.afrac(:);	
		fields = fieldnames(stats);
		for k=1:length(fields)
			if isfield(v,fields{k})
				spp.(sp).(fields{k})(j) = v.(fields{k});
			end
		end		
	end
end
% spp.fin.f is a list of cases matching the criteria for C. fin;
% spp.fin.f1 is the single best match
% stageFracW1 goes with f1, stageFracW goes with f, etc.

save -append runs_disko_stats spp	