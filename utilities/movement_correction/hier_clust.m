function [stt, stp, flag, scl, mc_flag] = hier_clust(acorr, Fs, pixs, scl, stype, m)
% find boundaries of stable/nonstable sections
%   Jinghao Lu, 05/15/2017

    [pixh, pixw, nf] = size(m, 'reg');
    if nargin < 4 || isempty(scl)
        defpar = default_parameters;
        scl = defpar.mc_scl;
    end
 
    if nargin < 5 || isempty(stype)
        defpar = default_parameters;
        ttype = defpar.ttype;
        stype = parse_type(ttype);
    end
    
    %% divide into sections %%
    mc_flag = true;
    flag = 1;
%     threst1 = hist_gauss(acorr, 0.1);
%     threst = scl * pixs;
% %     thres1 = 2 * hist_gauss(acorr, 0.5) - hist_gauss(acorr, 0.99);
%     thres1 = 2 * hist_gauss(acorr, 0.5);
    threst1 = prctile(acorr, 99);
    threst2 = prctile(acorr, 1);
    threst = scl * pixs;
%     thres1 = 2 * hist_gauss(acorr, 0.5) - hist_gauss(acorr, 0.99);
%     thres1 = 2 * prctile(acorr, 50);
    thres1 = movement_thres(acorr);
    thres1 = min(mad(acorr) + median(acorr), thres1);
%     if threst1 < threst || threst2 > threst
%         thres = thres1;
%     else
%         thres = threst; %%% no more than scl (percentage) of the image size %%%
%     end
    if threst1 > threst
        thres = thres1;
    else
        thres = threst; %%% no more than scl (percentage) of the image size %%%
    end
    scl = thres / pixs; %%% update new scl %%%
    ids = acorr > thres;
    ids = ~ids; %%% small dilation with small threshold: 0.2s each side %%%
%     ids = ~(imdilate(ids, strel('disk', round(Fs / 5)))); %%% small dilation with small threshold: 0.2s each side %%%
    [l, n] = bwlabeln(ids);
    stt = zeros(n, 1);
    stp = zeros(n, 1);
    for i = 1: n
        stt(i) = find(l == i, 1);
        stp(i) = find(l == i, 1, 'last') + 1;
    end
    
    if length(stt) == 1
        mc_flag = false;
    end
    
%     %% adjust long sections for balance %%
%     alen = stp - stt + 1;
%     na = length(alen);
%     idx = find(alen > na);
%     for i = 1: length(idx)
%         nt = ceil(alen(idx(i)) / na);
%         nstp = alen(idx(i)) / nt;
%         tmp = stt(idx(i)): nstp: stp(idx(i));
%         tmp = tmp(:);
%         stt = [stt; round(tmp(2: end - 1))];
%         stp = [stp; round(tmp(2: end - 1)) - 1];
%     end
%     stt = sort(stt);
%     stp = sort(stp);
    
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
    
    %% further adjust sections for memory fitness %%
    alen = stp - stt + 1;
    nff = max(alen);
    nsize = pixh * pixw * nff * stype; %%% size of single %%%
    nbatch = batch_compute(nsize);
    ebatch = round(nff / nbatch);
    idx = find(alen > ebatch);
    for i = 1: length(idx)
        tmp = stt(idx(i)): ebatch: stp(idx(i));
        tmp = tmp(:);
        stt = [stt; round(tmp(2: end - 1))];
        stp = [stp; round(tmp(2: end - 1)) - 1];
    end
    stt = sort(stt);
    stp = sort(stp);   
    
    %% final adjust sections for restricted length %%
    alen = stp - stt + 1;
    sstep = 10;
    for i = 1: length(alen)
        if alen(i) > sstep
            nc = ceil(alen(i) / sstep);
            sc = round(linspace(stt(i), stp(i), nc + 1));
            stt = [stt; sc(2: end - 1)' + 1];
            stp = [stp; sc(2: end - 1)'];
        end
    end
    stt = sort(stt);
    stp = sort(stp);   
end