function zsection_pad_stack(seis,depth,aflag,xdist,allx,nid)
% CSECTION(SEIS,BEG,DT,AFLAG,XDIST) plots a color seismogram 
% section of seismograms in 2D array SEIS. BEG is begin time 
% of section, DT is sample interval and AFLAG determines amplitude 
% scaling. If AFLAG < 0 (the default setting if AFLAG not specified)
% each trace is scaled to its maximum amplitude, if  AFLAG = 0 
% each trace is scaled to the maximum amplitude of the
% entire section, and if AFLAG > 0 then each trace is scaled to
% that amplitude (note this option allows direct comparison of 
% plots produced by different calls to SECTION). XDIST is vector
% containing spatial locations. If specified, traces are scaled
% and plotted to represent spatial distribution.
 
nx=length(allx);

ny=size(seis,1);
nt=size(seis,2);
if nargin < 4
  aflag=-1;
end

% Pad matrix xy
xymat=zeros(nx+1,nt);

% Set colormap.
r1=[(0:31)/31,ones(1,32)];
g1=[(0:31)/31,(31:-1:0)/31];
b1=[ones(1,32),(31:-1:0)/31];
rwb=[r1',g1',b1'];
colormap(rwb);

hold on

countb=zeros(nx,1); % will count number of waveforms in bin
% Normalize traces as specified.
if aflag < 0
  for ix=1:nx
      ind=find(xdist==allx(ix)); 
      if ~isempty(ind)
          countb(ix) = sum(nid(ind));
          stack = sum(seis(ind,:),1)./length(ind);
          xymat(ix,:)=stack/max(abs(stack)+0.000000001);
    %       text(xdist(iy),beg-1,num2str(nid(iy)))
      end
  end
elseif aflag == 0
  for ix=1:nx
      ind=find(xdist==allx(ix)); 
      if ~isempty(ind)
          countb(ix) = sum(nid(ind));
          stack = sum(seis(ind,:),1)./length(ind);
          xymat(ix,:)=stack/max(max(abs(seis))+0.000000001);
    %       text(xdist(iy),beg-1,num2str(nid(iy)))
      end
  end
else
  for ix=1:nx
      ind=find(xdist==allx(ix)); 
      if ~isempty(ind)
          countb(ix) = sum(nid(ind));
          stack = sum(seis(ind,:),1)./length(ind);
          xymat(ix,:)=stack/aflag;
    %       text(xdist(iy),beg-1,num2str(nid(iy)))
      end
  end
end
pcolor([allx;360],depth,xymat');
shading flat
xlim([0,363])
ylim([min(depth),max(depth)])
axis('ij');
if aflag > 0
  caxis([-aflag,aflag]);
else 
    caxis([-1,1]);
end
ylab=ylabel('Depth [km]');
xlab=xlabel('Backazimuth');
sprintf("zsection_pad: min=%f, max=%f",min(xymat(:)),max(xymat(:)))
% set(gca,'FontName','Helvetica','FontSize',16,'Clipping','off','layer','top');
% set(xlab,'FontName','Helvetica','FontSize',16);
% set(ylab,'FontName','Helvetica','FontSize',16);