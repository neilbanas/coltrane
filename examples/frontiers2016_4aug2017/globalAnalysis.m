load runs_global_stats stats cases p1 forcing1

dtP = unique(cases.dtP(:));
T0 = unique(cases.T0(:));
u0 = unique(cases.u0(:));
ToT = unique(cases.Tb0_over_T0(:));

clear dteff
for i = 1:length(dtP)
	forc = coltraneForcing('simple-g','dtP',dtP(i));
	dteff(i) = sum(forc.P./(1+forc.P));
end


u0s = u0([6 3 1]); % representative u0 for Cfin, Cgla, Chyp

figure
%f=cases.Tb0_over_T0==0.4 & stats.Fpot>=1;
f=cases.Tb0_over_T0==0.4 & stats.Fpot>=1 & ismember(cases.u0,u0s);
scatter(jitter(cases.T0(f),0.1),stats.Sapot(f),20,cases.u0(f),'o','filled');
set(gca,'yscale','log');
colormap(warmrainbow);
hold on
% some data. black = adult weight, blue = structural weight
plot(0.5,[140 370 1600],'ko'); % Swalethorp et al 2011
plot(0.5,[140*(1-0.32) 370*(1-0.34) 1600*(1-0.56)],'bo'); % Swalethorp et al 
plot(10.1,120,'ko'); % Peterson et al. 1986
plot([-1.2 1.9],252+65.*[-1 1],'k'); % Campbell et al. 2014
T=4:12;
PL = 2.87-0.036.*T; % Wilson et al. 2015
W = exp(3.57.*log(PL)+1.48); % Runge et al. 2006
plot(T,W,'k');
plot([4 8 12],[233 221 183],'ko-'); % Campbell et al. 2001
plot([15 15],[73 115],'ro'); % C. hel structural and total (Rey-Rassat et al 2002)
xlabel('T0');
ylabel('Sa pot');
axis square;
ylim([20 2000]);
xlim([-2 16]);
print -dpdf raw_figures/global_sizeVsTemp

figure
names = {'C fin','C gla','C hyp'};
for i=1:length(names)
	f = find(abs(u0-u0s(i))==min(abs(u0-u0s(i))));
	subplot(length(names),2,(i-1)*2+1);
	z = stats.ngenpot(:,:,1,f);
	z(stats.Fpot(:,:,1,f)==0) = nan;
	contourf(dteff,T0,z',[0 0.55:9.5],'edgecolor','none');
	caxis([0 4]);
	hold on
%	contour(dteff,T0,z',[0.55 55],'k');
	xl = xlim;
	xlim([0 xl(2)]);
	axis square;
	if i==1, title('ngenpot','interpreter','none'); end
	ylabel([names{i} '(u0=' num2str(u0s(i)) ')']);
	colorbar
	
	subplot(length(names),2,(i-1)*2+2);
	z = stats.Fpot_per_gen(:,:,1,f);
	z(stats.Fpot(:,:,1,f)==0) = nan;
	contourf(dteff,T0,z','edgecolor','none');
	caxis([0 4]);
	hold on
	contour(dteff,T0,z',[1 1],'w');
	xl = xlim;
	xlim([0 xl(2)]);
	axis square;
	if i==1, title('Fpot_per_gen','interpreter','none'); end
	colorbar
	
%	subplot(length(names),3,(i-1)*3+3);
	z1 = stats.Fpot_at1yr(:,:,1,f);
	z1(isnan(z1))=0;
	z2 = stats.Fpot_at2yr(:,:,1,f);
	z2(isnan(z2))=0;
	z3 = stats.Fpot_per_gen(:,:,1,f);
	z3(isnan(z3))=0;
	contour(dteff,T0,z1',[1 1],'r');
	hold on
	contour(dteff,T0,z2',[1 1],'b');
	contour(dteff,T0,z3',[1 1],'k');
	xl = xlim;
	xlim([0 xl(2)]);
	axis square;
	if i==1, title('viable at 1yr, 2yr, or at all'); end
	colorbar
end
colormap(warmrainbow);
orient tall

subplot(3,2,1);
plot(127,0.5,'ko');

subplot(3,2,3);
plot(127,0.5,'ko'); % Disko Bay (Madsen et al 2001)
plot([185 265],[11 11],'k'); % Newport (Huyer et al 2007 + Ks=1-3)
plot(40,0,'kx'); % Sheba (Ashjian et al 2003)
for i=1:20
	ebs_tice(i) = i*10-10;
	forc = coltraneForcing('bering','tice',ebs_tice(i),'P0IA',16,'dtbloom',20);
	ebs_dteff(i) = sum(forc.P./(3+forc.P));
	ebs_T0(i) = mean(forc.T);
end
plot(ebs_dteff,ebs_T0,'k'); % Bering Sea

subplot(3,2,5);
plot(127,0.5,'ko');


print -dpdf raw_figures/global_bySpecies
