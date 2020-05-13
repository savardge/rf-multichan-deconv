figure(100);
tstart=0;
tfin=40.0;
itb=fix((tstart)/ndt+1);
ite=fix((tfin)/ndt);
time=[0:size(dumr,2)-1].*ndt;

for k=1:size(resfunc,1)

% Plot radial and transverse sections
    t = time(itb:ite);
    if sum(abs(bdumr(k,itb:ite))) > 0
        
        subplot(2,2,1)
        mm = max(abs(bdumr(k,itb:ite)));
        plot(t,bdumr(k,itb:ite)./mm);title(sprintf("radial for k=%i",k))
        ylim([-1 1])
        hline(0,'k')
        
        subplot(2,2,3);
        plot(t,resfunc_r(k,itb:ite));title(sprintf("resfunc r for k=%i, beta = %f",k,betar(k)))
        
        subplot(2,2,2)
        mm = max(abs(bdumr(k,itb:ite)));
        plot(t,bdumt(k,itb:ite)./mm);title(sprintf("transverse for k=%i",k))
        ylim([-1 1])
        hline(0,'k')
        
        subplot(2,2,4);
        plot(t,resfunc_t(k,itb:ite));title(sprintf("resfunc t for k=%i, beta = %f",k,betat(k)))

        pause
        clf;
    end
end