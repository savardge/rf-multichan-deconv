
addpath('/data2/home/genevieve/research/TOOLBOXES/subtightplot')
subplot = @(m,n,p) subtightplot (m, n, p, [0.09 0.09], [0.1 0.1],[0.1 0.1]);
prefix='_bin'

% Determine ordering for plotting. Want to plot things in bins of
% back azimuth with each bin moving ordered in terms of increasing
% epicentral distance.
[dum,isort]=sort(fix(sbaz/10)+0.0001./spslow);

% Determine delays for Ps conversion from appropriate depth.
tshift=0;
bdumr0=zeros(npb2,nfft);bdumt0=zeros(npb2,nfft);
for ipb2=1:npb2
    wtshift=-sqrt(-1)*omega*tshift;
    xdumr=urft0(ipb2,:).*exp(-wtshift);
    xdumt=uzft0(ipb2,:).*exp(-wtshift);
    bdumr0(ipb2,:)=real(ifft([xdumr,fliplr(conj(xdumr(2:nfft/2)))],nfft));
    bdumt0(ipb2,:)=real(ifft([xdumt,fliplr(conj(xdumt(2:nfft/2)))],nfft));
end

%% Plot sections
% Pds.
lp=0.5;
hp=0.02
cparm=0.06%-1; %0.16;
tstart=0;
tfin=30.0;
% Plot radial and transverse sections
bdumr=fbpfilt(bdumr0,ndt,hp,lp,2,1);
bdumt=fbpfilt(bdumt0,ndt,hp,lp,2,1);
itb=fix((tstart-tshift)/ndt+1);
ite=fix((tfin-tshift)/ndt);

vbaz=0.5*(bzband(:,1)+bzband(:,2));

figure(1);clf;
set(gcf,'Color','w','Position',[192 46 1211 725])
subplot(2,1,1)
csection(bdumr(isort,itb:ite),tstart,ndt,cparm);
% csection(bdumr(isort,itb:ite),tstart,ndt,cparm,sbaz(isort));
% csection_pad(bdumr(isort,itb:ite),tstart,ndt,cparm,sbaz(isort),vbaz,nid(isort));
grid on;box on;set(gca,'layer','top')
title([sta, ', filter [',num2str(hp),'-',num2str(lp),'], radial component']); colorbar
subplot(2,1,2)
csection(-bdumt(isort,itb:ite),tstart,ndt,cparm)
% csection(-bdumt(isort,itb:ite),tstart,ndt,cparm,sbaz(isort))
% csection_pad(-bdumt(isort,itb:ite),tstart,ndt,cparm,sbaz(isort),vbaz,nid(isort));
grid on;set(gca,'layer','top')
title([sta, ', filter [',num2str(hp),'-',num2str(lp),'], transverse component']); colorbar
orient tall
if savefig==1
    export_fig(gcf,fullfile('./figures',[sta,'_Pds',prefix]),'-nocrop')
    saveas(gcf,fullfile('./figures',[sta,'_Pds',prefix,'.fig']))
end


% Ppds.
lp=0.75;
hp=0.02;
cparm=0.06%-1; %0.16;
% cparm=0.0
tstart=0;
tfin=30.0;
% Plot radial and transverse sections
bdumr=fbpfilt(bdumr0,ndt,hp,lp,2,1);
bdumt=fbpfilt(-bdumt0,ndt,hp,lp,2,1);
itb=fix((tstart-tshift)/ndt+1);
ite=fix((tfin-tshift)/ndt);

figure(2); clf;
set(gcf,'Color','w','Position',[192 46 1211 725])
subplot(2,1,1)
csection(bdumr(isort,itb:ite),tstart,ndt,cparm);
% csection_pad(bdumr(isort,itb:ite),tstart,ndt,cparm,sbaz(isort),vbaz,nid(isort));
grid on;set(gca,'layer','top')
title([sta, ', filter [',num2str(hp),'-',num2str(lp),'], radial component']); colorbar
subplot(2,1,2)
csection(bdumt(isort,itb:ite),tstart,ndt,cparm);
% csection_pad(bdumt(isort,itb:ite),tstart,ndt,cparm,sbaz(isort),vbaz,nid(isort));
grid on;set(gca,'layer','top')
title([sta, ', filter [',num2str(hp),'-',num2str(lp),'], transverse component']); colorbar
orient tall
if savefig==1
    export_fig(gcf,fullfile('./figures',[sta,'_Ppds',prefix]),'-nocrop')
    saveas(gcf,fullfile('./figures',[sta,'_Ppds',prefix,'.fig']))
end



%% Plot event distribution.
hl=figure(3);set(hl,'Color','w','Position',[192 46 1211 725])
clf
subplot(1,2,1)
yyaxis left
plot(sbaz(isort),'b+')
xlim([0,length(sbaz(isort))])
ylim([0,360]);
xlab=xlabel('Bin #');
ylab=ylabel('Back Azimuth [degrees]');
set(gca,'FontName','Helvetica','FontSize',16,'Clipping','off');
set(xlab,'FontName','Helvetica','FontSize',16);
set(ylab,'FontName','Helvetica','FontSize',16);

yyaxis right
plot(spslow(isort),'ro')
grid on;set(gca,'layer','top')
hold off
title(sta);

subplot(1,2,2)
% evmap % this is a polar plot with pslow as radius and azimuth
cn = 10; %ceil(max(nid)); % Number Of Colors
cm = colormap(jet(cn));
polarscatter(deg2rad(sbaz(isort)),spslow, [], cm(min([fix(nid(isort))',ones(size(nid(isort)))'.*cn],[],2),:),'filled','MarkerEdgeColor','k')
pax=gca;
pax.ThetaDir='clockwise';
pax.ThetaZeroLocation='top';
rticklabels({'p = 0.02','p = 0.04','p = 0.06','p = 0.08'})
cb=colorbar;
set(cb,'Ticks',[1:cn]*(1/cn),'TickLabels',[cellstr(num2str([1:cn-1]'));cellstr(['>=',num2str(cn)])])
cb.Label.String='# events';
tit=title('Event/Bin Distribution');
set(gca,'FontName','Helvetica','FontSize',11,'Clipping','off');
set(tit,'FontName','Helvetica','FontSize',16);
orient tall
if savefig==1
    export_fig(gcf,fullfile('./figures',[sta,'_bindistro',prefix]),'-nocrop')
    saveas(gcf,fullfile('./figures',[sta,'_bindistro',prefix,'.fig']))
end

eventmap;
if savefig==1
    export_fig(gcf,fullfile('./figures',[sta,'_bindistro_map',prefix]),'-nocrop')
    saveas(gcf,fullfile('./figures',[sta,'_bindistro_map',prefix,'.fig']))
end