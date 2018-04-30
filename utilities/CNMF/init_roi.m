function [roi, sig] = init_roi(data, sig, nIter)
% initialize roi spatial
%   Jinghao Lu, 06/11/2016

    if nargin < 3
        nIter = 1; 
    end
    [d1, d2, d3] = size(data);
    data = reshape(data, d1 * d2, d3);
    
    for iter = 1: nIter
        roi = max(data * sig, 0) / norm(sig);
        an = norm(roi);
        if an > 0
            roi = roi / an;
        else
            disp('degenerate component!')
        end
        sig = data' * roi;
    end
    roi = reshape(roi, d1, d2);
    sig = sig(:)';
end