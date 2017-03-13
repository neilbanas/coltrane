filename = '/Users/neil/datasets/coltrane_runs/disko/';
Nyears = 3;
paramCases = {'Nyears',Nyears,...
			  'u0',0.004:0.00025:0.012,...
			  'Ks',1,...
			  'm0',[0.06 0.07 0.08 0.09]};
forcingCases = {'dT',0};


% run the runs
if 1
	[cases,stats,p1,forcing1] = ...
		coltraneEnsemble('disko',forcingCases,paramCases,filename);
	save runs_disko_phi_stats stats cases p1 forcing1
else
	load runs_disko_phi_stats stats cases p1 forcing1
end


figure
loglog(stats.Sapot,stats.Fpot_at1yr,'-',stats.Sapot,stats.Fpot_at2yr,'-');
hold on
plot(xlim,[1 1],'k:');
xlabel('Sa pot');
ylabel('F');
xlim([50 1000]);
ylim([0.1 10]);
markdate([97 244 690]);
axis square
legend('0.06','0.07','0.08','0.09')
print -dpdf raw_figures/disko_F_vs_Sa


% three species cases in detail
p=p1;
p.m0 = 0.08;
forc0 = coltraneForcing('disko');
t0 = p.tspawn;

p.u0 = 0.01;
vfin = coltraneModel(forc0,p);
p.u0 = 0.007;
vgla = coltraneModel(forc0,p);
p.u0 = 0.005;
vhyp = coltraneModel(forc0,p);



