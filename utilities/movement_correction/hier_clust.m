function [stt, stp, flag, scl] = hier_clust(acorr, Fs, pixs, scl)
% find boundaries of stable/nonstable sections
%   Jinghao Lu, 05/15/2017

    if nargin < 4 || isempty(scl)
        defpar = default_parameters;
        scl = defpar.mc_scl;
    end
    
    %% divide into sections %%
    flag = 1;
%     thres1 = hist_gauss(acorr);
%     thres = min(scl * pixs, min(mad(acorr) + median(acorr), thres1)); %%% no more than scl (percentage0) of the image size %%%
    thres1 = gmm_bg(acorr);
    thres = min(scl * pixs, thres1); %%% no more than scl (percentage0) of the image size %%%
    scl = thres / pixs; %%% update new scl %%%
    ids = acorr > thres;
    ids = ~(imdilate(ids, strel('disk', round(Fs / 5)))); %%% small dilation with small threshold: 0.2s each side %%%
    [l, n] = bwlabeln(ids);
    stt = zeros(n, 1);
    stp = zeros(n, 1);
    for i = 1: n
        stt(i) = find(l == i, 1) + 1;
        stp(i) = find(l == i, 1, 'last') + 1;
    end
    if any(stt == 2)
        stt(stt == 2) = 1;
    end
    
    %% adjust long sections for balance %%
    alen = stp - stt;
    na = length(alen);
    idx = find(alen > na);
    for i = 1: length(idx)
        nt = ceil((stp(idx(i)) - stt(idx(i))) / na);
        nstp = (stp(idx(i)) - stt(idx(i))) / nt;
        tmp = stt(idx(i)): nstp: stp(idx(i));
        tmp = tmp(:);
        stt = [stt; round(tmp(2: end - 1))];
        stp = [stp; round(tmp(2: end - 1)) - 1];
    end
    stt = sort(stt);
    stp = sort(stp);
    
    %% generate fake sections %%
    if isempty(stt)
        nf = length(acorr) + 1;
        ns = (nf ^ 2 / 2) .^ (1 / 3); %%% best section numbers %%%
        ms = round(nf / ns);
        rg = (1: ms: nf)';
        if rg(end) < nf
            rg = [rg; nf];
        end
        stt = rg(1: 2: end);
        stp = rg(2: 2: end);
        luse = min(length(stt), length(stp));
        stt = stt(1: luse);
        stp = stp(1: luse);
        flag = 0;
    end
end