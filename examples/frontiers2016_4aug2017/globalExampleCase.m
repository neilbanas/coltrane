p = coltraneParams('Nyears',10,'u0',0.007,'Ks',1,'m0',0.08);
forc = coltraneForcing('simple-g','T0',1,'dtP',48);

[v,p,forcing] = coltraneModel(forc,p);
v.phi(v.phi==0) = nan;
v.a(isnan(v.D)) = nan;
v.pot.F0(v.pot.F0==0) = nan;
t = forcing.t(:,1);


figure
orient tall
subplot 511
contourf(t,p.tspawn,v.D',0:0.1:1,'edgecolor','none');
axis([0 3*365 0 365]);
set(gca,'xtick',(0:0.5:3)*365);
ylabel('D');
pbaspect([3 1 1]);
colorbar

subplot 512
contourf(t,p.tspawn,v.S',0:30:300,'edgecolor','none');
axis([0 3*365 0 365]);
set(gca,'xtick',(0:0.5:3)*365);
ylabel('S');
pbaspect([3 1 1]);
colorbar

subplot 513
v.phi(isnan(v.phi) & ~isnan(v.S)) = 0;
contourf(t,p.tspawn,v.phi',0:30:3000,'edgecolor','none');
axis([0 3*365 0 365]);
set(gca,'xtick',(0:0.5:3)*365);
ylabel('phi');
pbaspect([3 1 1]);
colorbar

subplot 514
contourf(t,p.tspawn,v.lnN'./log(10),-7:0.5:0,'edgecolor','none');
axis([0 3*365 0 365]);
set(gca,'xtick',(0:0.5:3)*365);
ylabel('log10 N');
pbaspect([3 1 1]);
colorbar

subplot 515
v.a(isnan(v.a) & ~isnan(v.S)) = 0;
contourf(t,p.tspawn,v.a',[0 0.5],'edgecolor','none');
axis([0 3*365 0 365]);
set(gca,'xtick',(0:0.5:3)*365);
ylabel('a');
pbaspect([3 1 1]);
colorbar
colormap(warmrainbow);

print -dpdf raw_figures/globalExample_statevars


figure
orient tall
subplot 511
plot(t,forcing.P(:,1),'k');
xlim([0 3*365]);
set(gca,'xtick',(0:0.5:3)*365);
ylabel('Prey');

subplot(5,1,[2 3]);
f = find(v.pot.Fopt > 0);
col = warmrainbow(2*length(f));
col = col(1:length(f),:);
for i=f
	coli = col(i-f(1)+1,:);
	plot(t,v.pot.F0(:,i),'color',coli);
	hold on
	plot(      p.tspawn(i).*[1 1],[0 4],'--','color',coli);
	plot(  365+p.tspawn(i).*[1 1],[0 4],'--','color',coli);
	plot(2*365+p.tspawn(i).*[1 1],[0 4],'--','color',coli);
	plot(  365+p.tspawn(i), v.pot.Fmult(1,i),'ko');
	plot(2*365+p.tspawn(i), v.pot.Fmult(2,i),'ko');
end
plot(xlim,[1 1],'k:');
axis([0 3*365 0 4]);
set(gca,'xtick',(0:0.5:3)*365);
ylabel('F');

subplot 514
for i=f
	coli = col(i-f(1)+1,:);
	plot(t,v.D(:,i),'color',coli);
	hold on
end
axis([0 3*365 0 1]);
set(gca,'xtick',(0:0.5:3)*365);
ylabel('D');

subplot 515
Fy = v.pot.Fopt.^(1/p.Nyears);
Fy(Fy==0) = nan;
plot(p.tspawn,Fy,'k');
ylim([0 2]);
xlim([0 3*365]);
set(gca,'xtick',(0:0.5:3)*365);
ylabel('Fopt per year');

print -dpdf raw_figures/globalExample_fitness
