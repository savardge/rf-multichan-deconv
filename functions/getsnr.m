function [snrall] = getsnr(fdir,zfiles,dtap,ndt0,ndec)

nfiles=length(zfiles);
snrall=zeros(nfiles,1);
for ifile=1:nfiles
    
    %% VERTICAL component
    
    kname=zfiles{ifile};
    fname=fullfile([fdir,'/OK/'],kname);   
    a=strsplit(kname,'.'); a=strjoin(a(5:9),'-');

   % Retrieve vertical component and header params.  
    head=readsac(fname);
    [tdum,dumz]=readsac(fname); tdum=tdum';dumz=dumz';
    beg=head.B;
    t1=head.T8;
    evdp=head.EVDP;
    % Define P-wave time window.
    % *** TO CHANGE WHEN TIME PICKS FIXED. meanwhile just use T7
    if (evdp > 100 && evdp <= 250 && head.T7>0) || isnan(head.T9) % then use "t7"
        t4=head.T7;
    else % then use "t9"
        t4=head.T9;
    end   
%     t4=min(hlength+beg-2*ndt0,t4);    
    
    it4(ifile)=round((t4-t1+dtap)/(ndt0*ndec));    
    its=round((t1-beg-dtap)/ndt0);
    itf=round((t4-beg)/ndt0);  
    ts=(t1-beg-dtap);
    tf=(t4-beg);  
    if isnan(its) || isnan(itf)
        kname
        t1
        t4
        head
        error('t1 or t4 incorrect.')
    end
           
    % Estimate noise on vertical component.
    dumz = detrendp(dumz, 0);
    dumz = detrendp(dumz, 1);
    dumz = bpfilt(dumz,ndt0,0.1,4);
    snrall(ifile) = std(dumz(its:itf))/std(dumz(1:its));
end
end

