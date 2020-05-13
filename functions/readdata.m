%----------------------------------------------------------
%
% SUBROUTINE READDATA
%
% Subroutine to open data files.
%
%-----------------------------------------------------


%% GET FILE LIST

% vertical
[~,out]=unix(['ls ',fdir,'/OK/',net,'.',sta,'*HZ.*SAC | xargs -i basename {}  ']);
zfiles=strsplit(out); zfiles=zfiles(1:end-1);
nfiles=length(zfiles);
% radial
[~,out]=unix(['ls ',fdir,'/OK/',net,'.',sta,'*HR.*SAC | xargs -i basename {}  ']);
rfiles=strsplit(out); rfiles=rfiles(1:end-1);
nfilesr=length(rfiles);
% transverse
[~,out]=unix(['ls ',fdir,'/OK/',net,'.',sta,'*HT.*SAC | xargs -i basename {}  ']);
tfiles=strsplit(out); tfiles=tfiles(1:end-1);
nfilest=length(tfiles);

if nfiles ~= nfilesr || nfiles ~= nfilest
    error('Non equal number of files.')
end

%% Preliminary calculations

% calculate differential times
% (e2p,tpp,pM,tM,p220,t220,p360,t360,p410,t410,p520,t520,p660,t660,ppmp,tdM,pcp,etc...)
difftimes;

% Choose length of FFT.
nfft=2^ceil(log(hlength/(ndecimate*ndt0))/log(2));
omega=(0:nfft/2)*2*pi/(nfft*ndt); % frequency array

%% INITIALIZATIONs

% Initialize arrays with zeros for faster input.
pft=zeros(nfiles,nfft); % P Fourier transform
tft=zeros(nfiles,nfft); % Transverse Fourier transform
sft=zeros(nfiles,nfft); % S Fourier transform
rft=zeros(nfiles,nfft); % radial Fourier transform
zft=zeros(nfiles,nfft); % Z Fourier transform

pwf=zeros(nfiles,nfft); % P waveform
swf=zeros(nfiles,nfft); % S waveform
twf=zeros(nfiles,nfft); % Transverse waveform

evlist=cell(nfiles,1); % filename list

evdp=zeros(1,nfiles); % event depth
evdist=zeros(1,nfiles); % event epidistance
pslow=zeros(1,nfiles); % p-wave slowness
baz=zeros(1,nfiles); % backazimuth
noisevert=zeros(nfiles,1); % noise on vertical component

%% Loop over file list and read
nover100under250=0; % counter for number of events between 100 and 250 km

for ifile=1:nfiles
    
     
    %% VERTICAL component
    
    kname=zfiles{ifile};
    fname=fullfile([fdir,'/OK/'],kname);   
    a=strsplit(kname,'.');a=strjoin(a(5:9),'-');
    evlist{ifile}=a;
    
   % Retrieve vertical component and header params.  
    head=readsac(fname);
    if ifile==1
        stalat=head.STLA;
        stalon=head.STLO;
        staelev=head.STEL;
    end
    [tdum,dumz]=readsac(fname); tdum=tdum';dumz=dumz';
    beg=head.B;
    evdp(ifile)=head.EVDP;
    evdist(ifile)=head.DIST;
    t1=head.T8;
    % Extract slowness (user0) and back azimuth for binning.
    pslow(ifile)=srad2skm(head.USER0);
    baz(ifile)=head.BAZ;
    
    % Define P-wave time window.
    % *** TO CHANGE WHEN TIME PICKS FIXED. meanwhile just use T9
    if (evdp(ifile) > 100 && evdp(ifile) <= 250 && head.T7>0) % then use "t4"
        if include100==0 % exclude events between 100-250 km
            continue
        end
%         t4=head.T9;
        t4=head.T7;
%         disp('100-250 km');
        nover100under250=nover100under250+1;

    else % then use "t3"
        t4=head.T9;
    end   
    t4=min(hlength+beg-2*ndt0,t4);    
    
    it4(ifile)=round((t4-t1+dtap)/(ndt0*ndec));    
    its=round((t1-beg-dtap)/ndt0);
    itf=round((t4-beg)/ndt0);  
    ts=(t1-beg-dtap);
    tf=(t4-beg);  
    if isnan(its) || isnan(itf)
        kname
        error('t1 or t4 incorrect.')
    end
    irpts=min(irpts,length(dumz)-itf); % Conflict with earlier definition???
            
    % Estimate noise on vertical component.
    bj=max(std(bpfilt(dumz(1:its),ndt0,0.1,1.0)),0.2);
        
    
    %decimate & detrend vertical channel trace
    
    zcomp=decimate(dumz(its:itf+irpts),ndec);
%     zcomp = detrend(zcomp);
%     zcomp=bpfilt(dumz(its:itf+irpts),ndt0,0.01,8);
%     zcomp=decimate(zcomp,ndec);
    zcomp=detrendp(zcomp, 0); % detrend (mean)
    zcomp=detrendp(zcomp, 1); % detrend (line fit)
%     zcomp=detrendp(zcomp, 2); % detrend (poly 2 fit)
    
    %% Radial component
    
    kname=rfiles{ifile};
    fname=fullfile([fdir,'/OK/'],kname);   
    
    % Retrieve radial component, note that we compensate here for
    % the fact that radial component is negative.        
%     rcomp=zeros(1,n2);    
    [tdum,dumr]=readsac(fname);  tdum=tdum';dumr=dumr';
    
    %decimate & detrend radial channel trace  
    rcomp=decimate(dumr(its:itf+irpts),ndec);
%     rcomp = detrend(rcomp);
%     rcomp=bpfilt(dumr(its:itf+irpts),ndt0,0.01,8);
%     rcomp=decimate(rcomp,ndec);
    rcomp=detrendp(rcomp, 0);   
    rcomp=detrendp(rcomp, 1);   
%     rcomp=detrendp(rcomp, 2);   
    
    %% Convert from R/Z to P/S on basis of slowness.
     % Pad vertical component
    zcomp=[zcomp,zeros(1,size(rcomp,2)-size(zcomp,2))];
    
    p2=pslow(ifile)*pslow(ifile);
    qa=sqrt(1/a02-p2);
    qb=sqrt(1/b02-p2);
    vpz=-(1-2*b02*p2)/(2*a0*qa);
    vpr=pslow(ifile)*b02/a0;
    vsr=(1-2*b02*p2)/(2*b0*qb);
    vsz=pslow(ifile)*b0;
    pcomp=vpz*zcomp-vpr*rcomp;
    scomp=vsz*zcomp-vsr*rcomp;
    % taper 
    pcomp=taper(pcomp,dtap,ndt,dtap,t4-t1+dtap);
    
    %% FOR TRANSVERSE COMPONENT (divide by 2 to remove free surface effect).
    
    kname=tfiles{ifile};
    fname=fullfile([fdir,'/OK/'],kname);      
    
    % Read transverse component
%     tcomp=zeros(1,n2);
    [tdum,dumr]=readsac(fname); tdum=tdum';dumr=dumr';
    
    %decimate & detrend radial channel trace   
    
    tcomp=decimate(dumr(its:itf+irpts),ndec);
%     tcomp = detrend(tcomp);
%     tcomp=bpfilt(dumr(its:itf+irpts),ndt0,0.01,8);
%     tcomp=decimate(tcomp,ndec);
    tcomp=detrendp(tcomp, 0);
    tcomp=detrendp(tcomp, 1);
%     tcomp=detrendp(tcomp, 2);
    tcomp = tcomp/2;
    
    %% Taper S and T waveforms (GS)
    tcomp=taper(tcomp,dtap,ndt,dtap,(length(tcomp)-1)*ndt-dtap); %GS
    scomp=taper(scomp,dtap,ndt,dtap,(length(scomp)-1)*ndt-dtap); %GS
    
    %% After rotation now truncate source (P).
    pcomp=[pcomp(1:it4(ifile)),zeros(1,nfft-it4(ifile))];
       
    % Compute fft and store in fft trace stack
    pft(ifile,:)=fft([pcomp,zeros(1,nfft-max(size(pcomp)))]/bj,nfft);
    sft(ifile,:)=fft([scomp,zeros(1,nfft-max(size(scomp)))]/bj,nfft);
    tft(ifile,:)=fft([tcomp,zeros(1,nfft-max(size(tcomp)))]/bj,nfft);
    rft(ifile,:)=fft([rcomp,zeros(1,nfft-max(size(rcomp)))]/bj,nfft);
    zft(ifile,:)=fft([zcomp,zeros(1,nfft-max(size(zcomp)))]/bj,nfft);
    % store waveforms:
    pwf(ifile,1:length(pcomp))=pcomp;
    swf(ifile,1:length(scomp))=scomp;
    twf(ifile,1:length(tcomp))=tcomp;
    noisevert(ifile)=bj;
    %% plot
%     figure(100);
%     subplot(3,1,1); section([pcomp],0,ndt,-1,['P'])
%     subplot(3,1,2); section([scomp],0,ndt,-1,['R'])
%     subplot(3,1,3); section([tcomp],0,ndt,-1,['T'])
%     
%     figure(101);        
%     [ppp,fp] = periodogram([pcomp,zeros(1,nfft-max(size(pcomp)))]/bj,[],[],100/ndec);
%     [pps,fs] = periodogram([scomp,zeros(1,nfft-max(size(scomp)))]/bj,[],[],100/ndec);
%     [ppt,ft] = periodogram([tcomp,zeros(1,nfft-max(size(tcomp)))]/bj,[],[],100/ndec);    
%     subplot(3,1,1); plot(fp,ppp); xlim([0 5])
%     subplot(3,1,2); plot(fs,pps); xlim([0 5])
%     subplot(3,1,3); plot(ft,ppt); xlim([0 5])
%     
%     figure(102)
%     plot(pcomp(1:length(scomp)))
%     hold on
%     plot(scomp,'r')
%     plot(tcomp,'g')
%     legend('P','S','T')
%     hold off
%     
%     pause
    
    
end
disp([num2str(nfiles), ' files read.'])
disp([num2str(nover100under250), ' events between 100-250 km.'])

% Get station coordinates.
stla=head.STLA;
stlo=head.STLO;




