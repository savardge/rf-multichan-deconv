% Examine bin
ii=16;
ibin=isort(ii);
ind=find(ipb==ipbuniq(ibin)); % find traces in bins
 
tstart=2.5;
tfin=10.0;

sbaz(ibin)
spslow(ibin)
        
% figure(100);
% subplot(1,2,1);section(pft(ind,1:nfft/2+1),0,ndt,-1); xlim([0,15]); title('FFT P wave')
% subplot(1,2,2);section(sft(ind,1:nfft/2+1),0,ndt,-1); xlim([0,15]); title('FFT S wave')

figure(101);
subplot(1,2,1);section(pwf(ind,:),0,ndt,-1); xlim([0,90]); title(' P wave')
subplot(1,2,2);section(swf(ind,:),0,ndt,-1); xlim([0,90]); title(' S wave')

figure(201)
subplot(2,1,1);
section(bdumr(ii,:),0,ndt,-1)
subplot(2,1,2);
section(bdumt(ii,:),0,ndt,-1)

csection(bdumt(ii,itb:ite),tstart,ndt,cparm);
