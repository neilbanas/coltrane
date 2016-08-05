load runs_disko_stats
stats.gl = stats.genlength./365;

% identify valid cases
f = stats.LEPn>1; 								% viable
redundant = zeros(size(cases.tegg));
for i=1:length(unique(cases.tegg(:)))-1
	redundant(i,:) = stats.LEPn(i,:)==stats.LEPn(i+1,:) & ...
					 stats.Wa(i,:)==stats.Wa(i+1,:);
end
f1 = find(~redundant & f & round(stats.gl)==1);	% f, broken down by gen length
f2 = find(~redundant & f & round(stats.gl)==2);
f3 = find(~redundant & f & round(stats.gl)==3);


% primary scatterplot of genlength vs Wa, colored by LEPn
figure
scatter(stats.Wa(f0),stats.gl(f0),50,[.7 .7 .7]);
hold on
scatter(stats.Wa(f),stats.gl(f),50,stats.LEPn(f),'filled');
colorbar;
caxis([1 6]);
colormap(tailedrainbow(0,0,3,6,0.01));
set(gca,'xscale','log');
xlim([50 3000]);
set(gca,'xtick',[50 100 200 500 1000 2000]);
ylim([0 4.5]);
xlabel('Wa (ugC)');
ylabel('gen length (years)');
title('LEPn');
markdate([spp.fin.W spp.gla.W spp.hyp.W]);
print -dpdf raw_figures/disko_genlength_Wa_LEPn


% same figure without colorscale, just valid/invalid LEPn
figure
plot(stats.Wa(f0),stats.gl(f0),'o','markeredgecolor',[.7 .7 .7]);
hold on;
plot(stats.Wa(f),stats.gl(f),'o');
set(gca,'xscale','log');
xlim([50 3000]);
set(gca,'xtick',[50 100 200 500 1000 2000]);
ylim([0 4.5]);
xlabel('Wa (ugC)');
ylabel('gen length (years)');
markdate([spp.fin.W spp.gla.W spp.hyp.W]);
print -dpdf raw_figures/disko_genlength_Wa


% plot annual cycle for single best case for each sp.	
for i=1:length(spp.names)
	sp = spp.names{i};
	fsp = spp.(sp).f1;
	figure
	subplot 311
	plot(p1.tspawn,cumsum(spp.(sp).stageFracW1([1:end 1],:)')');
	axis([0 365 0 1]);
	ylabel('Fraction of biomass');
	datetick('x','mmm','keeplimits');
	title(['best ' sp ': u0 = ' num2str(cases.u0(fsp)) ...
		   ', tegg = ' num2str(cases.tegg(fsp)) ...
		   ' (run ' num2str(fsp) ')']);
	subplot 312
	plot(p1.tspawn,spp.(sp).n1);
	ylabel('Egg production');
	xlim([0 365]);
	datetick('x','mmm','keeplimits');
	subplot 313
	plot(p1.tspawn,spp.(sp).afrac1([1:end 1]));
	axis([0 365 0 1]);
	ylabel('Active fraction');
	datetick('x','mmm','keeplimits');
	print('-dpdf',['raw_figures/disko_anncycle_best_' sp]);
end


% plot annual cycle for composite of all matching cases for each sp.
for i=1:length(spp.names)
	sp = spp.names{i};
	figure
	subplot 311
	stageFracWAvg = mean(spp.(sp).stageFracW,3);
	plot(p1.tspawn,cumsum(stageFracWAvg([1:end 1],:)')');
	axis([0 365 0 1]);
	ylabel('Fraction of biomass');
	datetick('x','mmm','keeplimits');
	title(['composite ' sp]);
	subplot 312
	nAvg = mean(spp.(sp).n,2);
	plot(p1.tspawn(:),nAvg);
	ylabel('Egg production');
	xlim([0 365]);
	datetick('x','mmm','keeplimits');
	subplot 313
	afracAvg = mean(spp.(sp).afrac,2);
	plot(p1.tspawn(:),afracAvg([1:end 1]));
	axis([0 365 0 1]);
	ylabel('Active fraction');
	datetick('x','mmm','keeplimits');
	print('-dpdf',['raw_figures/disko_anncycle_composite_' sp]);
end


% plot forcing
t = forcing1.t(:,1);
nn = 1:find(t==365);
figure
subplot 311
plot(t(nn),forcing1.Tb(nn,1),t(nn),forcing1.T(nn,1));
ylabel('Temperature');
xlim([0 365]);
datetick('x','mmm','keeplimits');
subplot 312
plot(t(nn),forcing1.P(nn,1));
ylabel('Prey');
xlim([0 365]);
datetick('x','mmm','keeplimits');
subplot 313
plot(t(nn),forcing1.P(nn,1)./(p1.Ks+forcing1.P(nn,1)));
ylabel('Prey saturation');
axis([0 365 0 1]);
datetick('x','mmm','keeplimits');
print -dpdf raw_figures/disko_anncycle_forcing


% Wa vs u0
figure
subplot 221
semilogx(stats.Wa(f1), cases.u0(f1),'o',...
		 stats.Wa(f2), cases.u0(f2),'o',...
		 stats.Wa(f3), cases.u0(f3),'o');
xlim([50 3000]);
set(gca,'xtick',[50 100 200 500 1000 2000]);
ylim([.004 .011]);
ylabel('u0');
xlabel('Wa');
axis square
print -dpdf raw_figures/disko_u0_vs_Wa

% various metrics vs Wa, colored by gen length
stats.tE50_jittered = jitter(stats.tE50,1);
stats.Ddiamin_jittered = jitter(stats.Ddiamin,0.01);
metrics = {'tE50_jittered','Ddiamin_jittered','capfrac','rhomean'};
figure
for i=1:length(metrics)
	subplot(2,2,i);
	semilogx(stats.Wa(f1),stats.(metrics{i})(f1),'o',...
			 stats.Wa(f2),stats.(metrics{i})(f2),'o',...
		 	 stats.Wa(f3),stats.(metrics{i})(f3),'o');
	xlim([50 3000]);
	set(gca,'xtick',[50 100 200 500 1000 2000]);
	xlabel('Wa');
	ylabel(metrics{i});
	axis square
end
% fix ranges
subplot(2,2,strmatch('capfrac',metrics));
ylim([0 1]);
subplot(2,2,strmatch('rhomean',metrics));
ylim([0 0.7]);
% add Swalethorp et al data (WE as fraction of total carbon) to rhomean
subplot(2,2,strmatch('rhomean',metrics));
hold on
plot(spp.fin.W.*[1 1],[0.32 0.32],'k*');
plot(spp.gla.W.*[1 1],[0.34 0.34],'k*');
plot(spp.hyp.W.*[1 1],[0.56 0.56],'k*');
print -dpdf raw_figures/disko_metrics_vs_Wa
% add stages to Ddiamin
subplot(2,2,strmatch('Ddiamin_jittered',metrics));
hold on
plot(xlim,(stage2D('C4')+stage2D('C5'))/2.*[1 1],'k--');
print -dpdf raw_figures/disko_variousMetrics