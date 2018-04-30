function [sn,psdx,ff] = get_noise_fft(Y,options)
% fft of all pixel traces, simplified from CNMF
%   Jinghao Lu, 06/11/2016
        
    defoptions = CNMFSetParms;
    if nargin < 2 || isempty(options); options = defoptions; end

    if ~isfield(options,'noise_range'); options.noise_range = defoptions.noise_range; end
    range_ff = options.noise_range;
    if ~isfield(options,'noise_method'); options.noise_method = defoptions.noise_method; end
    method = options.noise_method;

    dims = ndims(Y);
    sizY = size(Y);
    N = sizY(end);
    Fs = 1;        
    ff = 0: Fs / N: Fs / 2;
    indf=ff > range_ff(1);
    indf(ff > range_ff(2)) = 0;
    
    d = prod(sizY(1: dims - 1));
    Y = reshape(Y, d, N);    
    
    xdft = fft(Y, [], 2);
    xdft = xdft(:, 1: floor(N / 2) + 1);
    psdx = (1 / (Fs * N)) * abs(xdft) .^ 2;
    psdx(:, 2: end - 1) = 2 * psdx(:, 2: end - 1);
    switch method
        case 'mean'
            sn = sqrt(mean(psdx(:, indf) / 2, 2));
        case 'median'
            sn = sqrt(median(psdx(:, indf) / 2), 2);
        case 'logmexp'
            sn = sqrt(exp(mean(log(psdx(:, indf) / 2), 2)));
    end
end
