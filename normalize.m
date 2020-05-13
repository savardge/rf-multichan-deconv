% Here we perform normalization of the time series, and collect the
% normalized time series in bins.

% Filter parameters for std normalization.
nipbuniq=length(ipbuniq); % # of unique bins in p x baz matrix
sbaz=zeros(1,nipbuniq);
spslow=zeros(1,nipbuniq);
gcarc=zeros(1,nipbuniq);

% Initialization.
psumf1=zeros(1,nfft);
psumf2=zeros(1,nfft);
urft0=zeros(nbinp,nfft/2+1);
uzft0=zeros(nbinp,nfft/2+1);

pmax_r=zeros(nbinp,1);
pmax_t=zeros(nbinp,1);
pmax_pslow=zeros(nbinp,1);

dumr=zeros(nipbuniq,nfft);
dumz=zeros(nipbuniq,nfft);
wtshift=-sqrt(-1)*omega*tshift;
ipb2=0;

% Cycle through all bins.
nid = zeros(nipbuniq,1);
for ii=1:nipbuniq
    
    ind=find(ipb==ipbuniq(ii)); % find traces with same p and baz bins
   
    % If bin not empty then process.
    if ~isempty(ind)
        
        ipb2=ipb2+1;
        nid(ipb2)=length(ind); % number of traces stacked
                
        % Find mean p and baz for these bins
        [ip,ib]=ind2sub([nbinp,nbinbz],ipbuniq(ii));
        sbaz(ipb2)=mean(bzband(ib,:));
        spslow(ipb2)=mean(pband(ip,:));       
        isp(ipb2)=ip; % p for this RF
        isb(ipb2)=ib; % baz for this RF
        epid=interp1(e2p(:,2),e2p(:,1),spslow(ipb2));
        gcarc(ipb2) = km2deg(epid);
%         sprintf("num_traces= %02d, BAZ=%03d, p=%f, gcarc=%02d", nid(ipb2),sbaz(ipb2),spslow(ipb2),gcarc(ipb2))
        
        %% Simultaneous deconvolution
        % Use simdecf to produce deconvolved time series. Note that you
        % can send positive frequencies only and save time - the calculated
        % regularization parameter will be the same.
        %
        % UR/UZ
        %    [urft0(ipb2,:),psumf,beta]=simdecf(pft(ind,1:nfft/2+1),rft(ind,1:nfft/2+1),-5);
        %    [uzft0(ipb2,:)]=simdecf(pft(ind,1:nfft/2+1),zft(ind,1:nfft/2+1),beta);
        
        
        % SV/SH Note here that urft0 corresponds to radial and uzft0
        % corresponds to transverse - a bit confusing I admit.
        [urft0(ipb2,:),psumf1,betar(ipb2)]=simdecf(pft(ind,1:nfft/2+1),sft(ind,1:nfft/2+1),-5);
        [uzft0(ipb2,:),psumf2,betat(ipb2)]=simdecf(pft(ind,1:nfft/2+1),tft(ind,1:nfft/2+1),-5);
        %     [uzft0(ipb2,:)]=simdecf(pft(ind,1:nfft/2+1),tft(ind,1:nfft/2+1),beta);
        if betar(ipb2) == inf || betat(ipb2) == inf
            disp("No minimum found, skipping bin.")
            urft0(ipb2,:) = zeros(size(urft0(ipb2,:)));
            uzft0(ipb2,:) = zeros(size(uzft0(ipb2,:)));
%             continue
        end
        % Check initial data.
        rdat=real(ifft([rft(ind,1:nfft/2+1),fliplr(conj(rft(ind,2:nfft/2)))],nfft,2));
        zdat=real(ifft([zft(ind,1:nfft/2+1),fliplr(conj(zft(ind,2:nfft/2)))],nfft,2));                
        pdat=real(ifft([pft(ind,1:nfft/2+1),fliplr(conj(pft(ind,2:nfft/2)))],nfft,2));
        svdat=real(ifft([sft(ind,1:nfft/2+1),fliplr(conj(sft(ind,2:nfft/2)))],nfft,2));                
       
        % Determine resolution function amplitude for weighting. We do
        % this to correct for the effect of regularization on relative
        % amplitude.
        dump_r=real(ifft([psumf1,fliplr(conj(psumf1(2:nfft/2)))],nfft));
        pmax_r(ipb2)=max(dump_r);
        dump_t=real(ifft([psumf2,fliplr(conj(psumf2(2:nfft/2)))],nfft));
        pmax_t(ipb2)=max(dump_t);
%         sprintf("pmax_r = %f\t pmax_t = %f",pmax_r(ipb2),pmax_t(ipb2))
        resfunc_r(ipb2,:)=dump_r; 
        resfunc_t(ipb2,:)=dump_t;
        
        % Determine slowness dependent weighting (this empirical weighting
        % is based on numerical experiments using the ratio of Tu^(SP)/Tu^(PP),
        % that is the ratio of converted versus isomode transmission coefficients.
        % Set pmax=1 if you wish to maintain artificial TuSP/TuPP amplitude.
        pmax_pslow(ipb2)=(13.1*pslow(ipb2)-0.048);
        
        % Now determine time series
        urft0(ipb2,:)=urft0(ipb2,:).*exp(wtshift);
        uzft0(ipb2,:)=uzft0(ipb2,:).*exp(wtshift);
        dumr(ipb2,:)=...
            real(ifft([urft0(ipb2,:),fliplr(conj(urft0(ipb2,2:nfft/2)))],nfft));
        dumz(ipb2,:)=...
            real(ifft([uzft0(ipb2,:),fliplr(conj(uzft0(ipb2,2:nfft/2)))],nfft));


        %      xi=max(-dumr(ipb2,1:400))/max(-dumz(ipb2,1:400))
        %      svel(ipb2)=sqrt(1-1/sqrt(1+xi*xi))/(sqrt(2)*spslow(ipb2))
        
        %% plot
        
%         disp(['baz = ',num2str(sbaz(ipb2))])
%         disp(['p = ',num2str(spslow(ipb2))])
%         disp(['noise on zcomp = ',num2str(noisevert(ipb2))])
%         disp(['pmax radial = ',num2str(pmax_r(ipb2))])
%         disp(['pmax transv = ',num2str(pmax_t(ipb2))])
%         disp(['pmax slowness = ',num2str(pmax_pslow(ipb2))])
%         disp(['epidistance = ',num2str(epid)])
%         
%         figure(1);clf;
%         subplot(2,1,1);
%         section(rdat,0,ndt,-1); title('Radial')
%         subplot(2,1,2)
%         section(zdat,0,ndt,-1); title('Transverse')
%         
%         figure(2);clf;
%         subplot(3,1,1); section(pwf(ind,:),0,ndt,-1,['P']);
%         subplot(3,1,2); section(swf(ind,:),0,ndt,-1,['R']);
%         subplot(3,1,3); section(twf(ind,:),0,ndt,-1,['T']);
%         
%         figure(3);clf;
%         r_recfunc=fbpfilt(dumr(ipb2,:),ndt,0.03,0.75,2,1);
%         t_recfunc=fbpfilt(dumz(ipb2,:),ndt,0.03,0.75,2,1);
%         
%         subplot(2,1,1); plot([0:length(r_recfunc)-1].*ndt,r_recfunc); title('Radial RF')
%         xlim([0 30]); hold on; hline(0,'k-')
%         subplot(2,1,2); plot([0:length(t_recfunc)-1].*ndt,t_recfunc); title('Transverse RF')
%         xlim([0 30]); hold on; hline(0,'k-')
%         
%         pause
        
    end
end

npb2=ipb2;

% Determine ordering for plotting. Want to plot things in bins of
% back azimuth with each bin moving ordered in terms of increasing
% epicentral distance.
% [dum,isort]=sort(fix(sbaz/5)+0.0001./spslow);
[dum,isort]=sort(fix(sbaz)+spslow);
bazmid = sbaz(isort);
pmid = spslow(isort);
[bazmid_uniq,ic_baz]=unique(fix(bazmid));

% Determine delays for Ps conversion from appropriate depth.
tshift=0;
bdumr0=zeros(npb2,nfft);
bdumt0=zeros(npb2,nfft);
for ipb2=1:npb2
    wtshift=-sqrt(-1)*omega*tshift;
    xdumr=urft0(ipb2,:).*exp(-wtshift);
    xdumt=uzft0(ipb2,:).*exp(-wtshift);
    
    bdumr0(ipb2,:)=real(ifft([xdumr,fliplr(conj(xdumr(2:nfft/2)))],nfft));
%     bdumr0(ipb2,:) = bdumr0(ipb2,:)./pmax_r(ipb2);
%     bdumr0(ipb2,:)=detrend(bdumr0(ipb2,:),0);
%     bdumr0(ipb2,:)=detrend(bdumr0(ipb2,:),1);
    
    bdumt0(ipb2,:)=real(ifft([xdumt,fliplr(conj(xdumt(2:nfft/2)))],nfft));
%     bdumt0(ipb2,:) = bdumt0(ipb2,:)./pmax_t(ipb2);
%     bdumt0(ipb2,:)=detrend(bdumt0(ipb2,:),0);
%     bdumt0(ipb2,:)=detrend(bdumt0(ipb2,:),1);
end


