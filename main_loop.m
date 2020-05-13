stations=['SOLO';'TBO ';'GTOO';'SILO';'TIMO';'PNPO';'SUNO';'MATQ';'PKLO';...
    'KILO';'MALO';'WEMQ';'VLDQ';'VIMO';'EPLO';'CHGQ';'ULM ';'KAPO';'MUMO';'KASO'];
networks=[repmat('CN',7,1);'PO';'CN';'CN';'CN';'PO';repmat('CN',3,1);'PO';'CN';'CN';'CN';'CN'];
nsta = size(stations,1);

for is=20:nsta
    
    net = networks(is,:);
    sta = strtrim(stations(is,:));
    sprintf("Station %s.%s",net,sta)
    [stack,stats,p_window,svwf,shwf] = main2(net,sta);
    
    savefile = sprintf("RFdata_%s_%s.mat",net,sta);
    save(savefile, 'stack','stats','p_window','svwf','shwf')
    
end


[stack,stats,p_window,svwf,shwf] = main2('CN','MUMO');