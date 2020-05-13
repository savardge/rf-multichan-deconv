%% Make the grid.
% makegrid0; 
% . . Division in absolute slowness.
pmin=min(pslow) - 0.001;
pmax=max(pslow) + 0.001;
% avenumtr=20; % average number of trace per bin
nbinp=40; %round((1/avenumtr) * nfiles);
% dpstep=(pmax-pmin)/nbinp;
% pband=[(pmin:dpstep:pmax-dpstep);...
%     (pmin+dpstep:dpstep:pmax)]';
[idx,C]=kmeans(sort(pslow)',nbinp-1); C=sort(C);
countp=histc(idx,unique(idx));
disp(['Min # traces per p bin: ',num2str(min(countp))])
disp(['Max # traces per p bin: ',num2str(max(countp))])
disp(['Mean # traces per p bin: ',num2str(mean(countp))])
pband=[[pmin;C(1:end)],[C(1:end);pmax]];

% . . Division in back azimuth.
bzmin=0.0;
bzmax=360.0;
% bzstep=5.0;
% bzband=[(bzmin:bzstep:bzmax-bzstep);...
%     (bzmin+bzstep:bzstep:bzmax)]';
% nbinbz=size(bzband,1);
nbinbz=60;
[idx,C]=kmeans(sort(baz)',nbinbz-1); 
countb=histc(idx,unique(idx));
[C,ic]=sort(C); countb=countb(ic);
disp(['Min # traces per baz bin: ',num2str(min(countb))])
disp(['Max # traces per baz bin: ',num2str(max(countb))])
disp(['Mean # traces per baz bin: ',num2str(mean(countb))])
bzband=[[bzmin;C(1:end)],[C(1:end);bzmax]];


% . . Define slowness/baz grid.
npbz=nbinbz*nbinp;
bpbinmat=zeros(nfiles,nbinp,nbinbz);

for ifile=1:nfiles
%% Collect in appropriate bins, by storing bin ID of file. 
    %First find radial slowness index.
    ip=min(find(pslow(ifile) < pband(:,2))); % p index in pband array
    ib=min(find(baz(ifile) < bzband(:,2))); % baz index in bzband array
    ipb(ifile)=sub2ind([nbinp,nbinbz],ip,ib); % index in pbin x bzbin matrix
    bpbinmat(ifile,ip,ib)=1;
end

% Store number of unique grid indices in pbin x bzbin matrix
ipbuniq=unique(ipb);