if 1
	forcingCases = {'tice',0:5:190,...
						'T0mean',0:0.25:8,...
						'P0IA',16,...
						'dtbloom',20};
	paramCases = {'Nyears',2,...
				  'tegg',0,...
				  'u0',0.007,...
				  'Ks',3,...
				  'm0',0.08};
	filename = '~/datasets/coltrane_runs/bering/';

	[cases,stats,p1,forcing1] = ...
		coltraneEnsemble('bering',forcingCases,paramCases,filename);
	save runs_bering_stats stats cases p1 forcing1
else
	load runs_bering_stats stats cases p1 forcing1
end

if 1
	figure
	contourf(cases.tice,cases.T0mean,stats.Fpot_at1yr,0:0.2:3,...
		'edgecolor','none');
	hold on
	contour(cases.tice,cases.T0mean,stats.Fpot_at1yr,[1 1],'k');
	colormap(tailedrainbow(0,0.5,3,3));
	caxis([0 3]);
	colorbar;
	plot([75 75],[0 8],'k:');
	contour(cases.tice,cases.T0mean,stats.Fpot_per_yr - stats.Fpot_at1yr,[0.1 0.1],'w');
	moor = load('../../data/bering/BESTMAS_mooring_stats.mat');
	moor.tice(isnan(moor.tice)) = 0;
	pr = find(moor.years<2015);
	plot(moor.tice(pr,1),moor.T0mean(pr,1),'ro');
	plot(moor.tice(pr,4),moor.T0mean(pr,4),'bo');

	print -dpdf raw_figures/bering_fitness
end


if 1
	resol = 3;
	% M2, averages for 2003-05 and 2007-10
	forcingCases = {'tice',[0 105],...
					'T0mean',[4.5 6.2],...
					'P0IA',logspace(log10(8),log10(100),resol),...
					'dtbloom',20};
	paramCases = {'Nyears',2,...
				  'tegg',0,...
				  'u0',linspace(0.006,0.008,resol),...
				  'Ks',linspace(1,3,resol),...
				  'm0',logspace(log10(0.04),log10(0.16),resol)};
	filename = '~/datasets/coltrane_runs/bering_sens/';
	[cases,stats,p1,forcing1] = ...
		coltraneEnsemble('bering',forcingCases,paramCases,filename);
	save runs_bering_sens_stats stats cases p1 forcing1
else
	load runs_bering_sens_stats stats cases p1 forcing1
end
