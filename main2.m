% addpath('./functions')
% addpath(genpath('/Users/genevieve/TOOLBOXES/m_map'))

function [stack,stats,p_window,svwf,shwf] = main2(net,sta)

%     sta='SILO';
%     net='CN';
    fdir=['/Users/genevieve/metal_earth/MGBmethod/data/',sta]
    dec_fac = 5;
    ndt = 0.01 * dec_fac; % sampling after decimation
    
    % Read data
    disp('Reading data...')
    tic;
    [stats,p_window,~,svwf,shwf] = read_data(fdir,net,sta, dec_fac);
    toc
    
    % Define bins
    disp('Creating bins...')
    tic;
    bins = make_grid(stats);
    toc
    
    % Create new struct for RFs:
    stack = bins;
    clear bins
    
    % Delete bins that are empty
    stack = stack(~cellfun(@isempty,{stack.events_idx}));
    
    % Simultaneous time deconvolution
    disp("Simultaneous time deconvolution...")
    tic;
   
    for ib = 1:length(stack)
        
        % Get events in bin
        nevents = stack(ib).nevents;
        if nevents == 0; error("No events"); end
        id_evs = stack(ib).events_idx;
        
        % Deconvolve and get receiver function
        p_win = p_window(id_evs,:);
        sv_waveform = svwf(id_evs,:);
        sh_waveform = shwf(id_evs,:);
        bj = stats(id_evs).noise_vert;
        pslow = stack(ib).p;
        rf_sv = deconvolve(p_win,sv_waveform, bj, pslow, ndt);
        rf_sh = deconvolve(p_win,sh_waveform, bj, pslow, ndt);
        
        % assign to bins struct
        stack(ib).rf_sv = rf_sv;
        stack(ib).rf_sh = rf_sh;
        
        % Add a flag if no minimum was found in GCV
        if isempty(rf_sv); stack(ib).isBetaInf_sv = true; end
        if isempty(rf_sh); stack(ib).isBetaInf_sh = true; end

    end
    toc
    
    % plot

end

%% MAIN FUNCTIONS

function rf = deconvolve(pcomp, scomp, norm_fac, pslow, ndt)
% Simultaneous deconvolution
% Use simdecf to produce deconvolved time series. Note that you
% can send positive frequencies only and save time - the calculated
% regularization parameter will be the same.
        
    % Length for fft
    N = 200/ndt;
    nfft = 2^(nextpow2(N)+1);
    omega=(0:nfft/2)*2*pi/(nfft*ndt); % frequency array
    
    % Normalize
    if length(norm_fac) > 1; keyboard;end
    pcomp = diag(1./norm_fac) * pcomp;
    scomp = diag(1./norm_fac) * scomp;
    
    % Transform to frequency domain
    pft = fft(pcomp, nfft, 2);
    sft = fft(scomp, nfft, 2);
    
    % Apply simultaneous deconvolution
    [rf_ft,psum,beta]=simdecf(pft(:,1:nfft/2+1),sft(:,1:nfft/2+1),-1);
    
    if isinf(beta)
        % No minimum found, discard
        rf = [];
        return
    end
    
    % Scaling
    scale = 'none';
    switch scale
        case 'beta'
            % Determine resolution function amplitude for weighting. We do
        % this to correct for the effect of regularization on relative
        % amplitude.
            dump = real(ifft([psum,fliplr(conj(psum(2:nfft/2)))],nfft));
            pmax = max(dump);
        case 'pslow'
            % Determine slowness dependent weighting (this empirical weighting
        % is based on numerical experiments using the ratio of Tu^(SP)/Tu^(PP),
        % that is the ratio of converted versus isomode transmission coefficients.
        % Set pmax=1 if you wish to maintain artificial TuSP/TuPP amplitude.
            pmax = (13.1 * pslow - 0.048);
        case 'none'
            pmax = 1;
    end
    
    rf_ft = rf_ft * pmax;
  
    % Phase shift
    tshift = 0;
    wtshift = -sqrt(-1) * omega * tshift;
    rf_ft = rf_ft .* exp(-wtshift);
    
    % Transform back to time domain
    rf = real(ifft([rf_ft, fliplr(conj(rf_ft(2:nfft/2)))], nfft));

    if sum(rf) == 0
        error("RF is empty")
    end
end

function [stats,p_window,pwf,svwf,shwf] = read_data(fdir,net,sta, dec_fac)
    verbose = true;
    
    % Parameters
    ndt_raw = 0.01;
    taper_len = 3.0; 
    ndt = ndt_raw * 5;
    
    % Length of RF:
    rf_len_sec = 60;
    rf_len_samp = fix(rf_len_sec / ndt_raw);
    
    
    % Get file list
    zfiles = get_file_list(fdir, net, sta, 'Z');
    rfiles = get_file_list(fdir, net, sta, 'R'); 
    tfiles = get_file_list(fdir, net, sta, 'T'); 
    nevents = length(zfiles);
    
    % Initialize stats struct array
    stats = struct('Otime',[],...
                    'event_id',[],...
                    'baz',[],...
                    'pslow',[],...
                    'depth',[],...
                    'dist_km',[],...
                    'dist_deg',[],...
                    'magnitude',[],...
                    'fname_r',[],...
                    'fname_z',[],...
                    'fname_t',[],...
                    'tbegin',[],...
                    'p_win_start',[],...
                    'p_win_end',[],...
                    'p_predict',[],...
                    'SNR',[],...
                    'noise_vert',[],...
                    'is_z100_250',false,...
                    'include',false);
    stats = repmat(stats, nevents, 1);
    
    % Initialize waveform output
    p_window = zeros(nevents,rf_len_samp);
    pwf = zeros(nevents,rf_len_samp);
    shwf = zeros(nevents,rf_len_samp);
    svwf = zeros(nevents,rf_len_samp);

    % Loop over events and read data
    for iev = 1:nevents
        
        stats(iev).fname_z = zfiles{iev};
        stats(iev).fname_r = rfiles{iev};
        stats(iev).fname_t = tfiles{iev};
        
        % Get event info
        a=strsplit(stats(iev).fname_z,'.'); a = strjoin(a(5:9),'-');
        stats(iev).Otime = a;
        head = readsac(fullfile(fdir,'OK',stats(iev).fname_z));
        if round(head.DELTA,2) ~= 0.01
            disp(stats(iev).fname_z)
            sprintf("DELTA = %f", head.DELTA)
            error("Sampling of waveform is not 100 Hz!"); 
        end
        stats(iev) = extract_event_info(head, stats(iev));
        
        % Read in raw waveforms
        [tzraw,zraw] = readsac(fullfile(fdir,'OK',stats(iev).fname_z)); zraw = zraw';
        [~,rraw] = readsac(fullfile(fdir,'OK',stats(iev).fname_r)); rraw = rraw';
        [~,traw] = readsac(fullfile(fdir,'OK',stats(iev).fname_t)); traw = traw';
        N = min([length(zraw),length(traw),length(rraw)]);
        tzraw = tzraw(1:N); zraw = zraw(1:N); rraw = rraw(1:N); traw=traw(1:N);
        
        
        % Define window to use
        ts = stats(iev).p_win_start - stats(iev).tbegin - taper_len; % Start of P window - taper
        tf_pwin = stats(iev).p_win_end - stats(iev).tbegin; % end of P window
        its_raw = round(ts / ndt_raw);
        itf_raw = round(tf_pwin / ndt_raw);
        ite_raw = itf_raw + min(rf_len_samp, length(zraw) - itf_raw);  
        
        % Estimate noise on vertical component
        zraw_f = bpfilt(zraw,ndt_raw,0.1,1.0);
        zbefore = zraw_f(1:its_raw);
        zwin = zraw_f(its_raw:itf_raw);
        stats(iev).noise_vert =  max(std(bpfilt(zraw(1:its_raw),ndt_raw,0.1,1.0)),0.2);
        stats(iev).SNR = std(zwin)/std(zbefore);
        
        % Cut and pre-process
        dum = preprocess([zraw;rraw;traw], its_raw, ite_raw, dec_fac);
        zcomp = dum(1,:);
        rcomp = dum(2,:);
        tcomp = dum(3,:);        
        
        % Rotate RTZ to P-SV-SH
        [pcomp,svcomp,shcomp] = freetran(stats(iev).pslow, zcomp, rcomp, tcomp);
        
        % Taper
        % P window
        pwin_len_sec =  (tf_pwin - ts) + taper_len;
        pwin_len_samp = round(pwin_len_sec / ndt);
        pcomp_win = pcomp(1:pwin_len_samp); % cut
        tap_start = taper_len;
        tap_end = pwin_len_sec - taper_len - ndt; 
        pcomp_win = taper(pcomp_win,taper_len,ndt,tap_start,tap_end); % taper
 
        % S components
        tap_start = taper_len;
        tap_end = (length(svcomp)-1) * ndt - taper_len - ndt; 
        pcomp = taper(pcomp,taper_len,ndt,tap_start,tap_end);
        svcomp = taper(svcomp,taper_len,ndt,tap_start,tap_end);
        shcomp = taper(shcomp,taper_len,ndt,tap_start,tap_end);

        
        % plot
        if verbose == true
            figure(1);clf;
            % Raw
            t = tzraw;
            ax1 = subplot(3,2,1); plot(t,zraw);ylabel('Z','FontSize',18);title(sprintf('RAW : SNR = %f, noise vert = %f',stats(iev).SNR, stats(iev).noise_vert))
            ax3 = subplot(3,2,3); plot(t,rraw);ylabel('R','FontSize',18)
            ax5 = subplot(3,2,5); plot(t,traw);ylabel('T','FontSize',18)
            for ax = [ax1,ax3,ax5]
                axes(ax)
                v1 = vline([t(its_raw),t(ite_raw),t(itf_raw)],{'r','g','b'},{'start','RF','P_{window}'});
                v2 = vline([head.T0,head.T1,head.T5],'k',{head.KT0,head.KT1,head.KT5}); 
                set([v1,v2],'linewidth',2);
            end
            % Processed
            t = [0:length(pcomp)-1]*ndt;
            ax2 = subplot(3,2,2); plot(t,pcomp);ylabel('P','FontSize',18);title('Processed')
            hold on; plot(t(1:length(pcomp_win)),pcomp_win,'r','linewidth',2)
            ax4 = subplot(3,2,4); plot(t,svcomp);ylabel('SV','FontSize',18)
            ax6 = subplot(3,2,6); plot(t,shcomp);ylabel('SH','FontSize',18)
            pause;
        end
        
        % Assign
        p_window(iev,1:length(pcomp_win)) = pcomp_win; 
        pwf(iev,1:length(pcomp)) = pcomp;
        svwf(iev,1:length(pcomp)) = svcomp;
        shwf(iev,1:length(pcomp)) = shcomp;
        
        if sum(pcomp)==0 || sum(svcomp)==0 || sum(shcomp)==0
           keyboard 
        end
        
    end

end

function bins = make_grid(stats)
    
    % Slowness bins
    pmin=0.04;
    pmax=0.08;
    nbinp=50;
    dpstep=(pmax-pmin)/nbinp;
    pband=[(pmin:dpstep:pmax-dpstep);...
        (pmin+dpstep:dpstep:pmax)]';
    vp = mean(pband,2);
    
    % Backazimuth bins
    bzmin=0.0;
    bzmax=360.0;
    bzstep=5.0;
    bzband=[(bzmin:bzstep:bzmax-bzstep);...
        (bzmin+bzstep:bzstep:bzmax)]';
    vbaz = mean(bzband,2);
    nbinbaz = size(bzband,1);
    
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


    % Assembled bins in struct
    bins = struct('pmin',[],'pmax',[],'p',[],'bazmin',[],'bazmax',[],'baz',[],'events_idx',[],'nevents',0,'isempty',true);
    nbins = nbinp * nbinbaz;
    bins = repmat(bins, nbins, 1);
    
    for ib = 1:nbinbaz
        for ip = 1:nbinp
            
            % bin index
            ibin = sub2ind([nbinp,nbinbaz],ip,ib);
            
            % assign bin values
            bins(ibin).pmin = pband(ip,1);
            bins(ibin).pmax = pband(ip,2);
            bins(ibin).p = vp(ip);
            bins(ibin).bazmin = bzband(ib,1);
            bins(ibin).bazmax = bzband(ib,2);
            bins(ibin).baz = vbaz(ib);
            
            % find events in bin
            data_p = [stats.pslow];
            data_baz = [stats.baz];
            ind = find(data_baz >= bzband(ib,1) & data_baz <= bzband(ib,2)...
                & data_p >= pband(ip,1) & data_p <= pband(ip,2));
            
            % assign events
            if ~isempty(ind)
                bins(ibin).isempty = false;
                bins(ibin).events_idx = ind;
                bins(ibin).nevents = length(ind);
            else
                bins(ibin).isempty = true;
            end
            
        end
    end
    
end

%% HELPER FUNCTIONS

function traceOut = preprocess2(traceIn, ndt_raw, ts, te, taper_len, dec_fac)
    
    ntraces = size(traceIn,1);
    
    % Flip dimension to 1 trace per column
    traceIn = traceIn';
    
    % Window
    its = round((ts - taper_len) / ndt_raw);
    ite = min(round((te + taper_len) / ndt_raw), size(traceIn,1));
    
    % Cut
    traceIn = traceIn(its:ite,:);
 
    % Taper
    N = size(traceIn,1);
	win_percent = (2 * taper_len) / (N-1) * ndt_raw;
    trace = traceIn .* repmat(tukeywin(N, win_percent),1,ntraces); % 5% cosine window taper
    
    % Detrend
    trace = detrend(trace, 0); % demean
    trace = detrend(trace, 1); % detrend
    
    % Highpass filter
    fhigh = 0.05;
    trace = highpass(trace,fhigh,1/ndt_raw);
    
%     % decimate
    Nraw = size(trace,1);
    Ndec = ceil(Nraw/dec_fac);
    traceOut = zeros(Ndec,ntraces);
    for ii = 1:ntraces
        traceOut(:,ii) = decimate(trace(:,ii), dec_fac); % decimate
    end
    
    % Flip back to 1 trace per row
    traceOut = traceOut';
    
%     traceOut = fbpfilt(traceOut,ndt,0.05,4,4,1);
    
end

function traceOut = preprocess(traceIn, its, ite, dec_fac)
    
    ndt_raw = 0.01;
    ndt = ndt_raw * dec_fac;
    ntraces = size(traceIn,1);
    
    % Flip dimension to 1 trace per column
    traceIn = traceIn';
    
    % Cut
    trace = traceIn(its:ite,:);
 
    % Taper
%     N = size(traceOut,1);
% 	win_percent = 6 / (N-1)*ndt;
%     traceOut = traceOut .* repmat(tukeywin(N, win_percent),1,ntraces); % 5% cosine window taper
    
    % Detrend
    trace = detrend(trace, 0); % demean
    trace = detrend(trace, 1); % detrend
    
    % Highpass filter
    fhigh = 0.05;
    trace = highpass(trace,fhigh,1/ndt_raw);
    
    
%     % decimate
    Nraw = size(trace,1);
    Ndec = ceil(Nraw/dec_fac);
    traceOut = zeros(Ndec,ntraces);
    for ii = 1:ntraces
        traceOut(:,ii) = decimate(trace(:,ii), dec_fac); % decimate
    end
    
    % Flip back to 1 trace per row
    traceOut = traceOut';
    
%     traceOut = fbpfilt(traceOut,ndt,0.05,4,4,1);
    
end

function evout = extract_event_info(head, evin)
    evout = evin;
    evout.depth = head.EVDP;
    evout.dist_km = head.DIST;
    evout.dist_deg = head.GCARC;
    evout.pslow = srad2skm(head.USER0);
    evout.baz = head.BAZ;
    evout.magnitude = head.MAG;
    evout.p_predict = head.T0;
    
    % P window
    ts = head.T8;
    if (head.EVDP > 100 && head.EVDP <= 250 && head.T7>0) || isnan(head.T9)
        tf = head.T7;
        evout.is_z100_250 = true;
    else
        tf = head.T9;
        evout.is_z100_250 = false;
    end
    evout.p_win_start = ts;
    evout.p_win_end = tf;
    evout.tbegin = head.B;
end

function flist = get_file_list(fdir, net, sta, comp)
    tmp = dir(fullfile(fdir,'OK',[net,'.',sta,'*',comp,'.*SAC']));
    flist={tmp.name}';
end

function [pcomp,svcomp,shcomp] = freetran(pslow, zcomp, rcomp, tcomp)
    % Near surface crustal parameters for free surface correction.
    % Vp within 5.9-6.1, Vs within 3.7-3.8
    alpha = 6;
    beta = 3.75;
    
    % Transformation matrix parameters
    beta2 = beta^2; 
    alpha2 = alpha^2;
    p2 = pslow^2;
    qalpha = sqrt(1/alpha2 - p2);
    qbeta = sqrt(1/beta2 - p2);
    vpz = -(1 - 2 * beta2 * p2) / (2 * alpha * qalpha);
    vpr = pslow * beta2 / alpha;
    vsr = (1 - 2 * beta2 * p2) / (2 * beta * qbeta);
    vsz = pslow * beta;
    
    pcomp = vpz * zcomp - vpr * rcomp;
    svcomp = vsz * zcomp - vsr * rcomp;
    shcomp = tcomp ./ 2;
    
end

