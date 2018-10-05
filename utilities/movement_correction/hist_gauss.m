function [thres, wid, jud] = hist_gauss(sig, alf)
% Gauss fit of the 1d signal distribution
%   Jinghao Lu, 06/16/2016

    %% initialize %%
    if nargin < 2
        alf = 0.01;
    end
    nbins = round(length(sig) / 10);
    
    %% gauss fit %%
    [tmp1, ctrs] = hist(sig, nbins);
    tmp1 = tmp1(:);
%     [~, idm] = max(tmp1);
%     idm = max(3, idm);
    f = fit((1: nbins)', tmp1(:), 'gauss1');
%     f = fit((1: idm)', tmp1(1: idm), 'gauss1');
    wid = 2 * sqrt(log(1 / alf)) * f.c1;
    thres = wid / nbins * (max(sig) - min(sig)) + min(sig);
    
    if nargout > 2
        jud = [nbins / wid, r * nbins - f.b1];
    end
end