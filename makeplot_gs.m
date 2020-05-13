stations=['SOLO';'TBO ';'GTOO';'SILO';'TIMO';'PNPO';'SUNO';'MATQ';'PKLO';...
    'KILO';'MALO';'WEMQ';'VLDQ';'VIMO';'EPLO';'CHGQ';'ULM ';'KAPO';'MUMO';'KASO'];
networks=[repmat('CN',7,1);'PO';'CN';'CN';'CN';'PO';repmat('CN',3,1);'PO';'CN';'CN';'CN';'CN'];
nsta = size(stations,1);

for is=1:nsta


    fignum=is;

    ndt = 0.05;
    sta = strtrim(stations(is,:));
    net = networks(is,:);
    
    load(['RFdata_',net,'_', sta,'.mat' ])
    
    % Get RF info
    ix = find(~cellfun(@isempty,{stack.rf_sv}));
    stack2 = stack(ix);
    rf_sv = cell2mat({stack2.rf_sv}');
    rf_sh = cell2mat({stack2.rf_sh}');
    nevents = cell2mat({stack2.nevents}');
    mbaz = cell2mat({stack2.baz}');
    mpslow = cell2mat({stack2.p}');

    % filter
    lp = 0.05;
    hp = 1.2;
    rf_sv=fbpfilt(rf_sv,ndt,lp,hp,2,1);
    rf_sh=fbpfilt(-rf_sh,ndt,lp,hp,2,1);


    % cut
    itb=fix((0)/ndt+1);
    ite=fix((20)/ndt);
    rf_svc = rf_sv(:,itb:ite);
    rf_shc = rf_sh(:,itb:ite);

    cparm = -1;
    figure(fignum); clf;
    set(gcf,'Color','w','Position',[192 46 1211 725])
    subplot(2,1,1)
    csection2(rf_svc,0,ndt,cparm,nevents);
    grid on;set(gca,'layer','top')
    title([sta, ', SV component']); 
    % colorbar
    % set(gca,'XTick',ic_baz,'XTickLabel',bazmid_uniq,'TickDir','out')
    xtickangle(270)
    xlabel('Backazimuth')

    subplot(2,1,2)
    csection2(rf_shc,0,ndt,cparm,nevents);
    grid on;set(gca,'layer','top')
    % set(gca,'XTick',ic_baz,'XTickLabel',bazmid_uniq,'TickDir','out')
    xtickangle(270)
    xlabel('Backazimuth')
    title([sta, ', SH component']); 
    colorbar
end
%%
% 
for k=1:size(rf_svc,1)
    
    figure(1);clf;
    nev = stack2(k).nevents
    t = [0:size(svwf,2)-1]*ndt;
    for ie=1:nev
        id = stack2(k).events_idx(ie);
    	subplot(nev+1,1,ie)
        plot(t,p_window(id,:),'b',t,svwf(id,:),'r','linewidth',2)
    end
    a2=subplot(nev+1,1,nev+1);
    t = [0:size(rf_svc,2)-1]*ndt;
    plot(t,rf_svc(k,:),'r','linewidth',2)
    title(sprintf("p = %5.3f, baz = %03d, # events = %d",mpslow(k),mbaz(k),nevents(k)))
%     linkaxes([a1,a2],'x');
    pause;
end

%%
% 
% for k=1:100
%     sprintf("event ID = %i",stack2(k).events_idx)
%     [bj(k),bj2(k)]
%     t = [0:size(bdumr,2)-1]*ndt;
%     figure(1);clf;
%     a1=subplot(2,1,1);
% %     plot(pwf(k,:),'linewidth',2);hold on
% %     plot(p_window1(k,:)','linestyle','--','linewidth',2);
%     plot(t,bdumr(k,:),'b','linewidth',2)
%     title(sprintf("original: p = %5.3f, baz = %03d, # events = %d",spslow(k),sbaz(k),nid(k)))
%     a2=subplot(2,1,2);
% %     plot(swf(k,:),'linewidth',2);hold on
% %     plot(svwf1(k,:)','linewidth',2,'linestyle','--');
%     plot(t,rf_sv(k,:),'r','linewidth',2)
%     hold on;plot(t,bdumr(k,:),'b','linewidth',0.5)
%     title(sprintf("new: p = %5.3f, baz = %03d, # events = %d",mpslow(k),mbaz(k),nevents(k)))
%     linkaxes([a1,a2],'x');
%     pause;
% end
% 
% idx = 117
% figure(2);clf;
% subplot(2,1,1);
% plot(p_window1(idx,:))
% hold on
% plot(pwf(idx,:))
% subplot(2,1,2)
% plot(svwf1(idx,:))
% hold on
% plot(swf(idx,:))
% 
% idx = 117
% figure(2);clf;
% subplot(2,1,1);
% plot(real(fft(p_window1(idx,:),nfft)))
% hold on
% plot(real(pft(idx,:)))
% subplot(2,1,2)
% plot(real(fft(svwf1(idx,:),nfft)))
% hold on
% plot(real(sft(idx,:)))