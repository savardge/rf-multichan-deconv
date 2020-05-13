% Specify Azimuthal projection with a radius of 100 degrees from
% station.
figure(100);
clf
set(gcf,'Color','w')
maxed=100;

% Projection.
m_proj('Azimuthal Equidistant','lon',stlo,'lat',stla,...
    'rad',maxed,'rectbox','circle');

% Plot coastlines.
m_coast;

% Plot grid.
m_grid('box','on','xtick',3,'ytick',3);
hold on;

% Get corners of plot so that one can scale small circles.
xyaxis=axis;
rad=xyaxis(2);

% Plot station
[xsta,ysta]=m_ll2xy(stlo,stla);
plot(xsta,ysta,'k^','MarkerFaceColor','r')
text(xsta,ysta,sta,'FontSize',14)

% Calculate earthquake locations
for ipb2=1:nipbuniq
    
    % Determine backazimuth and epicentral distance
    % of node locations.
    %    sbaz=mean(bzband(isp(ipb2),isb(ipb2),:));
    %    spslow=mean(pband(isp(ipb2),:));
    sepd(ipb2)=interp1(e2p(:,2),e2p(:,1),spslow(ipb2));
    dd(ipb2)=km2deg(sepd(ipb2));
    [latout,lonout] = reckon(stla,stlo,dd(ipb2),sbaz(ipb2));
    [xpt,ypt]=m_ll2xy(lonout,latout);
    
    % Revert to plot coordinates.
%     phi=sbaz(ipb2)*pi/180;
%     xpt=rad*sepd(ipb2)*sin(phi)/maxed;
%     ypt=rad*sepd(ipb2)*cos(phi)/maxed;
    plot(xpt,ypt,'r.')
    hold on;
%     plot(xpt,ypt,'r+')
%     hold on;
    text(xpt,ypt,num2str(nid(ipb2)));
%     nid(ipb2);
end

%hold off
load station_coords.mat
[xpt,ypt]=m_ll2xy(metalearthstations.lon,metalearthstations.lat);
plot(xpt,ypt,'k.')

%%
ref1 = [52.294792, -164.409409] % tip of Alaskan panhandle
ref2 = [34.556507, 126.555180] % South Korea
m_line([ref1(2),stlo], [ref1(1),stla],'linestyle','-','linewidth',2,'color','r')
hold on
m_line([ref2(2),stlo], [ref2(1),stla],'linestyle','-','linewidth',2,'color','r')

% SE corridor
ref1 = [-36.569737, -53.624591] % east S.Am.
ref2 = [-4.618318, -88.074152] % west S.Am.
m_line([ref1(2),stlo], [ref1(1),stla],'linestyle','-','linewidth',2,'color','b')
hold on
m_line([ref2(2),stlo], [ref2(1),stla],'linestyle','-','linewidth',2,'color','b')

