%% Make the grid.
% makegrid0; 
% . . Division in absolute slowness.
% pmin=min(pslow) - 0.001;
% pmax=max(pslow) + 0.001;
pmin=0.04;
pmax=0.08;
% avenumtr=20; % average number of trace per bin
nbinp=50; %round((1/avenumtr) * nfiles);
dpstep=(pmax-pmin)/nbinp;
pband=[(pmin:dpstep:pmax-dpstep);...
    (pmin+dpstep:dpstep:pmax)]';

% . . Division in back azimuth.
bzmin=0.0;
bzmax=360.0;
bzstep=5.0;
bzband=[(bzmin:bzstep:bzmax-bzstep);...
    (bzmin+bzstep:bzstep:bzmax)]';

% bazv = 2:3:360;
% bzband = nan(length(2:3:360),2);
% for k=1:length(bazv)
%     [bazv(k)-2, bazv(k),bazv(k)+2]
%     bzband(k,:) = [bazv(k)-2, bazv(k)+2];
% end


% % % NW corridor
% ref1 = [52.294792, -164.409409]; % tip of Alaskan panhandle
% ref2 = [34.556507, 126.555180]; % South Korea
% bzmin = floor(azimuth(stla,stlo,ref1(1),ref1(2)));
% bzmax = ceil(azimuth(stla,stlo,ref2(1),ref2(2)));
% [bzmin, bzmax];
% bzband = [bzmin, bzmax; bzmax, bzmin];

% SE corridor
% ref1 = [-36.569737, -53.624591]; % east S.Am.
% ref2 = [-4.618318, -88.074152]; % west S.Am.
% bzmin = floor(azimuth(stla,stlo,ref1(1),ref1(2)));
% bzmax = ceil(azimuth(stla,stlo,ref2(1),ref2(2)));
% bzband = [bzmin, bzmax];

nbinbz=size(bzband,1);

% . . Define slowness/baz grid.
npbz=nbinbz*nbinp;
bpbinmat=zeros(nfiles,nbinp,nbinbz);

ipb=nan(nfiles,1);
for ifile=1:nfiles
%% Collect in appropriate bins, by storing bin ID of file. 
    %First find radial slowness index.
    ip=min(find(pslow(ifile) < pband(:,2))) % p index in pband array
%     ib=min(find(baz(ifile) < bzband(:,2))); % baz index in bzband array
    ib=find(baz(ifile) < bzband(:,2) & baz(ifile) >= bzband(:,1))
    if isempty(ib); continue; end
    ipb(ifile)=sub2ind([nbinp,nbinbz],ip,ib); % index in pbin x bzbin matrix
    bpbinmat(ifile,ip,ib)=1;
end

% Store number of unique grid indices in pbin x bzbin matrix
ipbuniq=unique(ipb(~isnan(ipb)));
