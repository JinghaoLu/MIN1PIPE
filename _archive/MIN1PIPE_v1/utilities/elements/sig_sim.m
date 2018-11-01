function sigcorfn = sig_sim(datasmthf, cutofff, pkcutofff, flag)
% calculate signal similarity matrix 
%   Jinghao Lu, 06/10/2016

    if nargin < 4
        flag = true;
    end
    
    nseed = size(datasmthf, 1);
    sigcor = zeros(nseed, nseed);
    sigcort = zeros(nseed, nseed);
    Fc = [eps, 0.5]; % Cut-off frequencies
    Fn = 10 / 2;
    Wn = Fc / Fn;
    [B, A] = butter(1, Wn, 'bandpass');
    datafilter = hilbert(filter(B, A, datasmthf'));
    aar = real(datafilter);
    aai = imag(datafilter);
    phaseall = (atan(aai ./ aar))';
    if flag %%% really compute the time domain correlation %%%
        for i = 1: nseed
            %%% find valid period for seeds correlation %%%
            pxluse = vld_prd_slct(datasmthf(i, :), cutofff(i), pkcutofff(i));
            
            %%% find useful time series of all the traces %%%
            datacorr = phaseall(i, pxluse);
            dataother = phaseall(:, pxluse);
            
            %%% compute 2 similarities %%%
            tp1 = datasmthf(i, pxluse);
            tp2 = datasmthf(:, pxluse);
            sigcort(i, :) = corr(tp1', tp2');
            sigcor(i, :) = norm_inner(datacorr, dataother');
        end
    else
        sigcor = corr(datasmthf');
    end
    tt = max(sigcor, sigcort);
    sigcorfn = min(tt, tt');
end