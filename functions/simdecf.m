function [rft,xft,betax] = simdecf(wft,vft,betan)

% SIMDECF Simultaneous deconvolution of multiple seismograms in
% the frequency domain. Inputs are wavelet estimates WFT, data
% VFT (both in frequency domain and of dimension M X N where
% M is number of seismograms and N is number of frequencies),
% and regularization parameter BETAN. If BETAN < 0 then
% an optimum parameter BETAX is sought using Generalized Cross
% Validation and used to produce impulse response RFT, and model
% resolution kernel XFT.

% Calculate dimensions of wft, and assume that minimum dimension is
% number of traces, maximum size is number of frequencies.
nm = min(size(wft));
nn = max(size(wft));

% Compute denominators.
if nm == 1
    wwft = wft.*conj(wft);
    vwft = vft.*conj(wft);
else
    wwft = sum(wft .* conj(wft));
    vwft = sum(vft .* conj(wft));
end

% If nonzero betan provided by user, use it for deconvolution.
% Otherwise compute best beta using Generalized Cross Validation.
if betan > 0
    betax = betan;
else
%     beta = exp(-1000 : 0.2 : 40);
%     beta = exp(-16:0.2:16);
    beta = exp(-900 : 0.2 : 40);
    gcvf = zeros(1, length(beta));
    for ib = 1:length(beta)
        
        % Define operator W W* / (W W* + B) and deconvolve to get impulse response in
        % frequency domain.
        wwft2 = wwft + beta(ib);
        rft = vwft ./ wwft2;
        xft = wwft ./ wwft2;
        
        % Compute model norm.
        %modnorm(ib)=norm(rft)^2;
        
        % Compute data misfit. Note misfit is numerator of GCV function.
        % Note also earlier mistake where norm(nft)^2 was norm(nft).
        if nm == 1
            nft = vft - wft .* rft;
            misfit = norm(nft)^2;
        else
            misfit = 0.0;
            for im = 1 : nm
                nft = vft(im,:) - wft(im,:) .* rft;
                misfit = misfit + norm(nft)^2;
            end
        end
        
        % Compute denominator and GCV function.
        den = nn * nm - real(sum(xft));
        den = den * den;
        gcvf(ib) = misfit / den;
    end
    
    % Compute best beta.
    [~, ibest] = min(gcvf);
    betax = beta(ibest);
%     figure(100);clf
%     plot(1:length(gcvf),gcvf,'r');hold on;plot(ibest,gcvf(ibest),'ko','linewidth',2)
%     title(sprintf("beta = %f, min = %f, max=%f",betax, exp(-16),exp(16)))
%     pause;
        
    % If minimum not found inform user.
    if ibest == 1 || ibest == length(beta)
        sprintf('WARNING: No minimum found for GCV')
        disp('index at minimum and no of seismograms');
        disp('change search limits')
%         [ibest,nm]
        betax = inf;
        rft = zeros(1,size(wft,2));
        xft = rft;
%         figure(100);clf
%         plot(beta,gcvf,'r');
%         error("No minimum found")
        return
    else
        % Final estimate.
        wwft2 = wwft + betax;
        rft=vwft ./ wwft2;
        xft=wwft ./ wwft2;
        return
    end
    
end

