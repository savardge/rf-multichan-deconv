%% Single-depth mapping

%% use the AK135 continental model
% [vmod.z, ~, vmod.vp, vmod.vs, ~, ~ ] =  ak135( 'cont' );

%% Shaeffer and Lebedev Vs model: SL2013NA
% Schaeffer, A. J. and S. Lebedev. Imaging the North American continent using waveform inversion of global and USArray data. 
% Earth and Planetary Science Letters, 402, pg 26-41, 2014. doi:10.1016/j.epsl.2014.06.014.

fid=fopen(['/Users/genevieve/metal_earth/velocity_model/SL2013NA_1D_0.25D/',sta,'.mod'],'r');
C = textscan(fid,'%f %f %f');
fclose(fid);
vmod.z=flipud(C{1});
vmod.vp=flipud(C{2});
vmod.vs=flipud(C{3});
ix=find(vmod.z>=0);
vmod.z=vmod.z(ix);
vmod.vp=vmod.vp(ix);
vmod.vs=vmod.vs(ix);

%% Plot radial and transverse sections
lp=0.8;
lp1=0.03;
bdumr = fbpfilt(bdumr0,ndt,lp1,lp,2,1);
bdumt = fbpfilt(bdumt0,ndt,lp1,lp,2,1);

tshift=0
tstart = 0;%2.5;
tfin = 40.0;
itb = fix((tstart-tshift)/ndt+1);
ite = fix((tfin-tshift)/ndt);
Prfr0 = bdumr(isort,itb:ite);
Prft0 = -bdumt(isort,itb:ite);

%% Collapse slowness bins: pick 1 per backazimuth bin
vbaz=0.5*(bzband(:,1)+bzband(:,2));

collapse = false
if collapse
    nx = length(vbaz);
    Prfr = zeros(nx+1, size(Prfr0,2));
    Prft = zeros(nx+1, size(Prft0,2));
    parr = zeros(nx+1, 1);
    bazarr = parr;

    xdist = sbaz(isort);
    allx = vbaz;
    nid = nid(isort);
    countb=zeros(nx,1);
    for ix=1:nx
        ind=find(xdist==allx(ix));
        if ~isempty(ind)
            [count,imax] = max(nid(ind));
            if count > 0
                countb(ix) = count;
                ikeep = ind(imax);
                Prfr(ix,:) = Prfr0(ikeep,:);
                Prft(ix,:) = Prft0(ikeep,:);
                parr(ix) = pmid(ikeep);
                bazarr(ix) = bazmid(ikeep);
            end
        end
    end
else
    Prfr = Prfr0;
    Prft = Prft0;
    parr = pmid;
    bazarr = bazmid;
    countb = nid;
end

%% convert stacked seismograms to depth

Prfr = num2cell( Prfr  ,2);
Prft = num2cell( Prft ,2);
time = [0:(length(Prfr{1})-1)] * ndt + tstart;
time = repmat(num2cell(time,2),size(Prfr,1),1);

dz = 0.25; % 0.05
maxz = 300;
maxxy = 150;

% [ depthr, eposr, nposr, seisr ] = mapPsSeis2depth_1d( time, Prfr, pmid, bazmid, dz, vmod.z, vmod.vp, vmod.vs );
% [ deptht, epost, npost, seist ] = mapPsSeis2depth_1d( time, Prft, pmid, bazmid, dz, vmod.z, vmod.vp, vmod.vs );
[ depthr, eposr, nposr, seisr ] = mapPsSeis2depth_1d( time, Prfr, parr, bazarr, dz, vmod.z, vmod.vp, vmod.vs );
[ deptht, epost, npost, seist ] = mapPsSeis2depth_1d( time, Prft, parr, bazarr, dz, vmod.z, vmod.vp, vmod.vs );


%% Now plot depth section
aflag = -1;

%RADIAL
[~,ncols] = cellfun(@size, seisr);
z=depthr{1}; z=z(1:min(ncols));
nt=min(ncols);
tmp=cellfun(@(x) x(1:nt)',seisr, 'UniformOutput',false);
PRFrmat=cell2mat(tmp)';

figure(figno1);clf; set(gcf,'Color','w')
subplot(2,2,1)
zsection(PRFrmat,z,aflag);
title(['Radial component, filtered between ',num2str(lp1),' - ',num2str(lp),' Hz'])
grid on;box on; 
set(gca,'layer','top')
set(gca,'XTick',ic_baz,'XTickLabel',bazmid_uniq,'TickDir','out')
xtickangle(270)
xlabel('Backazimuth')
ylim([0,250])
% set(gca,'Position',[0.1300    0.1100    0.7750    0.8150])
if savefig==1
    export_fig(gcf,fullfile('./figuresz',[sta,'_lp',replace(num2str(lp),'.','-'),'_radial_bin']),'-nocrop','-m2')
%     saveas(gcf,fullfile('./figuresz',[sta,'_lp',replace(num2str(lp),'.','-'),'_radial_baz']))
end

% figure(13);clf
subplot(2,2,3);
zsection_pad_stack(PRFrmat,z,aflag,bazarr,vbaz,countb)
title(['Radial component, filtered between ',num2str(lp1),' - ',num2str(lp),' Hz'])
grid on;box on;set(gca,'layer','top')
% set(gca,'Position',[0.1300    0.1100    0.7750    0.8150])
if savefig==1
    export_fig(gcf,fullfile('./figuresz',[sta,'_lp',replace(num2str(lp),'.','-'),'_radial_bin']),'-nocrop','-m2')
%     saveas(gcf,fullfile('./figuresz',[sta,'_lp',replace(num2str(lp),'.','-'),'_radial_baz']))
end
ylim([0,250])



%TRANSVERSE
[~,ncols] = cellfun(@size, seist);
z=deptht{1}; z=z(1:min(ncols));
nt=min(ncols);
tmp=cellfun(@(x) x(1:nt)',seist, 'UniformOutput',false);
PRFtmat=cell2mat(tmp)';

% figure(figno2);clf;
subplot(2,2,2)
zsection(PRFtmat,z,aflag);
title(['Transverse component, filtered between ',num2str(lp1),' - ',num2str(lp),' Hz'])
grid on;box on;set(gca,'layer','top')
set(gca,'XTick',ic_baz,'XTickLabel',bazmid_uniq,'TickDir','out')
xtickangle(270)
xlabel('Backazimuth')
ylim([0,250])
% set(gca,'Position',[0.1300    0.1100    0.7750    0.8150])
if savefig==1
    export_fig(gcf,fullfile('./figuresz',[sta,'_lp',replace(num2str(lp),'.','-'),'_transv_bin']),'-nocrop','-m2')
%     saveas(gcf,fullfile('./figuresz',[sta,'_lp',replace(num2str(lp),'.','-'),'_transv_bin']))
end

% figure(15);clf
subplot(2,2,4)
zsection_pad_stack(PRFtmat,z,aflag,bazarr,vbaz,countb)
title(['Transverse component, filtered between ',num2str(lp1),' - ',num2str(lp),' Hz'])
grid on;box on;set(gca,'layer','top')
% set(gca,'Position',[0.1300    0.1100    0.7750    0.8150])
if savefig==1
    export_fig(gcf,fullfile('./figuresz',[sta,'_lp',replace(num2str(lp),'.','-'),'_transv_baz']),'-nocrop','-m2')
%     saveas(gcf,fullfile('./figuresz',[sta,'_lp',replace(num2str(lp),'.','-'),'_transv_baz']))
end
ylim([0,250])



%% Cones
% figure(2); clf;
% plotRFgather_rayp( depth, seis , pmid, bazmid, numel(seis) , 1, true)
% ax = axis();
% axis( [ ax(1), ax(2), ax(3), maxz] );
% ylabel('Depth (km)')
% depth=depthr;
% seis=seisr;
% epos=eposr; npos=nposr;
% 
% % combine position data
% allpos = [ [epos{:}]', [npos{:}]', [depth{:}]' ];
% allv = [seis{:}]';
% 
% % restrict depth and horizontal
% idx2 = find( allpos(:,3)<maxz & ( allpos(:,1).^2 + allpos(:,2).^2 )<maxxy^2 );
% allpos = allpos( idx2, : );
% allv = allv(idx2);

% % plot in 3d
% figure(4); clf;
% maxA = 0.5*max(max(abs(allv)));
% clims = [-maxA, maxA ];
% mycmap = repmat( [linspace(0,1,32), linspace(1,0,32)]', 1, 3 );
% mycmap(1:32,3) = ones(32,1);
% mycmap(33:64,1) = ones(32,1);
% colormap(  mycmap );
% scatter3( allpos(:,1), allpos(:,2), allpos(:,3), 20*ones(numel(allv),1 ), allv, ...
%     'filled', 'LineWidth', 1 ); hold on;
% set(gca, 'CLim', clims );
% set(gca,'ZDir','reverse');
% xlabel( 'E-W distance (km)')
% ylabel( 'N-S distance (km)')
% zlabel( 'Depth (km)' );
% axis equal;

%% Cones

% RADIAL
% depth=depthr;
% seis=seisr;
% epos=eposr; npos=nposr;
% % combine position data
% allpos = [ [epos{:}]', [npos{:}]', [depth{:}]' ];
% allv = [seis{:}]';
% % restrict depth and horizontal
% idx2 = find( allpos(:,3)<maxz & ( allpos(:,1).^2 + allpos(:,2).^2 )<maxxy^2 );
% allpos = allpos( idx2, : );
% allv = allv(idx2);
% [lat,lon]=SDC2(allpos(:,1),allpos(:,2),1, stla, stlo);
% fid=fopen([sta,'_cone_radial.xyzv'],'w');
% fprintf(fid,'%f\t%f\t%f\t%f\t%f\t%f\n',[lat,lon,allpos(:,1),allpos(:,2),allpos(:,3),allv]');
% fclose(fid);

% figure(10); clf;
% maxA = 0.5*max(max(abs(allv)));
% clims = [-maxA, maxA ];
% mycmap = repmat( [linspace(0,1,32), linspace(1,0,32)]', 1, 3 );
% mycmap(1:32,3) = ones(32,1);
% mycmap(33:64,1) = ones(32,1);
% colormap(  mycmap );
% scatter3( allpos(:,1), allpos(:,2), allpos(:,3), 20*ones(numel(allv),1 ), allv, ...
%     'filled', 'LineWidth', 1 ); hold on;
% set(gca, 'CLim', clims );
% set(gca,'ZDir','reverse');
% xlabel( 'E-W distance (km)')
% ylabel( 'N-S distance (km)')
% zlabel( 'Depth (km)' );
% axis equal;


% TRANSVERSE
% depth=deptht;
% seis=seist;
% epos=epost; npos=npost;
% % combine position data
% allpos = [ [epos{:}]', [npos{:}]', [depth{:}]' ];
% allv = [seis{:}]';
% % restrict depth and horizontal
% idx2 = find( allpos(:,3)<maxz & ( allpos(:,1).^2 + allpos(:,2).^2 )<maxxy^2 );
% allpos = allpos( idx2, : );
% allv = allv(idx2);
% [lat,lon]=SDC2(allpos(:,1),allpos(:,2),1, stla, stlo);
% fid=fopen([sta,'_cone_transv.xyzv'],'w');
% fprintf(fid,'%f\t%f\t%f\t%f\t%f\t%f\n',[lat,lon,allpos(:,1),allpos(:,2),allpos(:,3),allv]');
% fclose(fid);

% figure(11); clf;
% maxA = 0.5*max(max(abs(allv)));
% clims = [-maxA, maxA ];
% mycmap = repmat( [linspace(0,1,32), linspace(1,0,32)]', 1, 3 );
% mycmap(1:32,3) = ones(32,1);
% mycmap(33:64,1) = ones(32,1);
% colormap(  mycmap );
% scatter3( allpos(:,1), allpos(:,2), allpos(:,3), 20*ones(numel(allv),1 ), allv, ...
%     'filled', 'LineWidth', 1 ); hold on;
% set(gca, 'CLim', clims );
% set(gca,'ZDir','reverse');
% xlabel( 'E-W distance (km)')
% ylabel( 'N-S distance (km)')
% zlabel( 'Depth (km)' );
% axis equal;