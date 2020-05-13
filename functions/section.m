function [h]=section(seis,beg,dt,aflag,stalst,ic)

% SECTION(SEIS,BEG,DT,AFLAG) plots a seismogram section of 
% seismograms in 2D array SEIS. BEG is begin time of section, DT
% is sample interval and AFLAG determines amplitude scaling. If
% AFLAG < 0 each trace is scaled to its maximum amplitude (default), 
% if AFLAG = 0 each trace is scaled to the maximum amplitude of the
% entire section, and if AFLAG > 0 then each trace is scaled to
% that amplitude (note this option allows direct comparison of 
% plots produced by different calls to SECTION).

% goodcolors
% black=[0 0 0];
% col=[black;darkred;red;orange;yellow;lightgreen;darkgreen;lightblue;darkblue;purple;pink];
col=colormap(lines(10));

ny=size(seis,1);
nt=size(seis,2);
if nargin < 4
  aflag=-1;
  stalst='';
  ic=1;
elseif nargin < 6
    ic=1;
end

yaxe=[1:ny];
if aflag < 0
  for iy=1:ny
    xymat(iy,:)=iy-0.7*seis(iy,:)/max(abs(seis(iy,:))+0.0000001);
  end
elseif aflag == 0
  for iy=1:ny
    xymat(iy,:)=iy-8.0*seis(iy,:)/max(max(abs(seis))+0.0000001);
  end
else
  for iy=1:ny
    xymat(iy,:)=iy-1.0*seis(iy,:)/aflag;
  end
end
time=[0:nt-1]*dt+beg;
%plot(time,xymat,[col(ic),'-']);
h=plot(time,xymat,'k');
set(h,'Color',col(ic,:))
axis('ij');
xlab=xlabel('Time [s]');
ylab=ylabel('Trace #');
ylim([0,ny+1]);
set(gca,'YTick',(1:ny))
if exist('stalst','var')
    set(gca,'YTickLabel',stalst)
end
set(gca,'FontName','Helvetica','FontSize',14,'Clipping','off','layer','top');
set(xlab,'FontName','Helvetica','FontSize',14);
set(ylab,'FontName','Helvetica','FontSize',14);

