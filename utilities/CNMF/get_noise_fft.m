function [sn,psdx,ff] = get_noise_fft(m,options)
% fft of all pixel traces, simplified from CNMF
%   Jinghao Lu, 06/11/2016
        
    defoptions = CNMFSetParms;
    if nargin < 2 || isempty(options); options = defoptions; end

    if ~isfield(options,'noise_range'); options.noise_range = defoptions.noise_range; end
    range_ff = options.noise_range;
    if ~isfield(options,'noise_method'); options.noise_method = defoptions.noise_method; end
    method = options.noise_method;

    [pixh, pixw, nf] = size(m, 'reg');
    Fs = 1;        
    ff = 0: Fs / nf: Fs / 2;
    indf = ff > range_ff(1);
    indf(ff > range_ff(2)) = 0;
    
    d = pixh * pixw;
    stype = parse_type(m.reg(1, 1, 1));
    nsize = d * nf * stype; %%% size of single %%%
    nbatch = batch_compute(nsize);
    ebatch = ceil(d / nbatch);
    eb = floor(ebatch / pixh);
    idbatch = [1: eb: pixw, pixw + 1];
    nbatch = length(idbatch) - 1;

    sn = zeros(d, 1);
    for i = 1: nbatch
        tmp = m.reg(1: pixh, idbatch(i): idbatch(i + 1) - 1, 1: nf);
        tmp = reshape(tmp, pixh * (idbatch(i + 1) - idbatch(i)), nf);
        xdft = fft(tmp, [], 2);
        xdft = xdft(:, 1: floor(nf / 2) + 1);
        psdx = (1 / (Fs * nf)) * abs(xdft) .^ 2;
        psdx(:, 2: end - 1) = 2 * psdx(:, 2: end - 1);
        switch method
            case 'mean'
                sn(pixh * (idbatch(i) - 1) + 1: pixh * (idbatch(i + 1) - 1)) = sqrt(mean(psdx(:, indf) / 2, 2));
            case 'median'
                sn(pixh * (idbatch(i) - 1) + 1: pixh * (idbatch(i + 1) - 1)) = sqrt(median(psdx(:, indf) / 2), 2);
            case 'logmexp'
                sn(pixh * (idbatch(i) - 1) + 1: pixh * (idbatch(i + 1) - 1)) = sqrt(exp(mean(log(psdx(:, indf) / 2), 2)));
        end
    end
end
