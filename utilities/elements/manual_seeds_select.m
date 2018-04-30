function [roi, sig, idusef, bg, bgf, datasmthf, cutofff, pkcutofff] = manual_seeds_select(frame, Fs, sz)
% manually select seeds to use, and then auto-initialize
%   Jinghao Lu, 02/06/2018
    
    %% initialization %%
    %%% initialize parameters %%%
    if nargin < 2 || isempty(Fs)
        defpar = default_parameters;
        Fs = defpar.Fsi_new;
    end    
    
    if nargin < 3 || isempty(sz)
        defpar = default_parameters;
        sz = defpar.neuron_size;
    end
    
    %% manual input of seeds %%
    [pixh, pixw, nf] = size(frame);
    figure(1)
    clf
    imagesc(max(frame, [], 3));
    [x, y] = ginput;
    
    %% data preparation %%
    idusef = sub2ind([pixh, pixw], round(y), round(x));
    datusef = reshape(frame, pixh * pixw, nf);
    datusef = datusef(idusef, :);
    gsw = 1 * Fs;
    datasmthf = convn(datusef, gausswin(gsw)' / sum(gausswin(gsw)), 'valid');
    st = sort(datasmthf, 2);
    bs = zeros(size(st));
    for i = 1: size(st, 1)
        bs(i, :) = linspace(st(i, 1), st(i, end), size(st, 2));
    end
    stcr = st - bs;
    [~, minp] = min(stcr, [], 2);
    [~, maxp] = max(stcr, [], 2);
    minids = sub2ind(size(st), (1:size(st,1))', minp);
    maxids = sub2ind(size(st), (1:size(st,1))', maxp);
    cutofff = (st(minids) + st(maxids)) / 2;
    pkcutofff = st(minids);

    %% initialize roi and sig %%
    res = frame;
    nseed = length(idusef);
    roi = zeros(pixh * pixw, nseed);
    sig = zeros(nseed, nf);
    swin = 3 * sz;
    cthres = 0.9; %%% previous 0.95 %%%
    
    for i = 1: nseed
        %%% get the current seed position %%%
        [x, y] = ind2sub([pixh, pixw], idusef(i));
        
        %%% prepare small patch of data around the seed %%%
        rg = [max(1, x - swin), min(pixh, x + swin); max(1, y - swin), min(pixw, y + swin)];
        lsml = prod(diff(rg, 1, 2) + 1);
        ctr = [x - rg(1, 1) + 1, y - rg(2, 1) + 1];
        a = res(rg(1, 1): rg(1, 2), rg(2, 1): rg(2, 2), :);
        sa = frame(rg(1, 1): rg(1, 2), rg(2, 1): rg(2, 2), :);
                
        %%% pixel correlation %%%
        tuse = vld_prd_slct(datasmthf(i, :), cutofff(i), pkcutofff(i));

        %%% coarse spatial init %%%
        saa = reshape(sa, lsml, nf);
        daa = saa(:, tuse);
        datt = datusef(i, :);
        dat = datt(tuse);
        ind = norm_inner(daa, dat') > cthres;
        ind = reshape(ind, diff(rg(1, :)) + 1, diff(rg(2, :)) + 1);
        [l, ~] = bwlabeln(ind);
        ii = l(ctr(1), ctr(2));
        msk = l == ii;
        data = bsxfun(@times, sa, msk);
        
        %%% get temporal and (spatial) init; previous method faster %%%
        [atmp, ctmp] = init_roi(data(:, :, tuse), dat');
        
        %%% rescale the atmp and ctmp %%%
        rt = max(dat) / max(ctmp);
        ctmp = datusef(i, :) / rt;
        
        %%% update rois and sigs %%%
        tmp  = zeros(pixh, pixw);
        tmp(rg(1, 1): rg(1, 2), rg(2, 1): rg(2, 2)) = atmp;
        roi(:, i) = reshape(tmp, pixh * pixw, 1);
        sig(i, :) = ctmp;
        tmp = saa - reshape(atmp, lsml, 1) * ctmp;
        
        %%% update remaining dataset and find next %%%
        res(rg(1, 1): rg(1, 2), rg(2, 1): rg(2, 2), :) = min(a, reshape(tmp, diff(rg(1, :)) + 1, diff(rg(2, :)) + 1, nf));
        if mod(i, round(nseed / 10)) == 0
            disp(['Done #', num2str(i), '/', num2str(nseed)])
        end
    end
    
    roi = sparse(roi);
      
    %%% get background roi and sig (unnecessary) %%%
    resmax = max(reshape(res, pixh * pixw, nf), 0);
    W0 = mean(resmax, 2);
    H0 = mean(resmax, 1);
    bg = W0;
    bgf = H0;
    disp('Done manual pix select')
end