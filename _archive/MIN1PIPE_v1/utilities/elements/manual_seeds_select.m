function [roi, sig, idusef, bg, bgf, datasmthf, cutofff, pkcutofff] = manual_seeds_select(m, Fs, sz)
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
    [pixh, pixw, nf] = size(m, 'reg');
    stype = parse_type(class(m.reg(1, 1, 1)));
    nsize = pixh * pixw * nf * stype; %%% size of single %%%
    nbatch = batch_compute(nsize);
    ebatch = ceil(nf / nbatch);
    idbatch = [1: ebatch: nf, nf + 1];
    nbatch = length(idbatch) - 1;
    imaxf = zeros(pixh, pixw);
    for i = 1: nbatch
        tmp = m.reg(1: pixh, 1: pixw, idbatch(i): idbatch(i + 1) - 1);
        imaxf = max(cat(3, max(tmp, [], 3), imaxf), [], 3);
    end
    figure(1)
    clf
    imagesc(imaxf);
    hold on
    key = 'NaN';
    x = [];
    y = [];
    while ~isempty(key)
        [xt, yt, key] = ginput(1);
        plot(xt, yt, '.r')
        x = [x, xt];
        y = [y, yt];
    end
    hold off
    close
    
    %% data preparation %%
    idusef = sub2ind([pixh, pixw], round(y), round(x));
    datusef = zeros(length(idusef), nf);
    stype = parse_type(class(m.reg(1, 1, 1)));
    nsize = pixh * pixw * nf * stype; %%% size of single %%%
    nbatch = batch_compute(nsize);
    ebatch = ceil(nf / nbatch);
    idbatch = [1: ebatch: nf, nf + 1];
    nbatch = length(idbatch) - 1;
    imeantf = zeros(1, nf);
    imeanf = zeros(pixh, pixw);
    for i = 1: nbatch
        tmp = m.reg(1: pixh, 1: pixw, idbatch(i): idbatch(i + 1) - 1);
        imeanf = (imeanf * (idbatch(i) - idbatch(1)) + double(sum(tmp, 3))) / (idbatch(i + 1) - idbatch(1));
        tmp = reshape(tmp, pixh * pixw, idbatch(i + 1) - idbatch(i));
        imeantf(idbatch(i): idbatch(i + 1) - 1) = double(mean(tmp, 1));
        datusef(:, idbatch(i): idbatch(i + 1) - 1) = tmp(idusef, :);
%         disp(num2str(i))
    end
    
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
    nseed = length(idusef);
    roi = zeros(pixh * pixw, nseed);
    sig = zeros(nseed, nf);
    swin = 2 * sz;
    cthres = 0.9; %%% a high enough value for pixel clustering as initialization, not very important %%%
        
    for i = 1: nseed
        %%% get the current seed position %%%
        [x, y] = ind2sub([pixh, pixw], idusef(i));
        
        %%% prepare small patch of data around the seed %%%
        rg = [max(1, x - swin), min(pixh, x + swin); max(1, y - swin), min(pixw, y + swin)];
        lsml = prod(diff(rg, 1, 2) + 1);
        ctr = [x - rg(1, 1) + 1, y - rg(2, 1) + 1];
        sa = m.reg(rg(1, 1): rg(1, 2), rg(2, 1): rg(2, 2), :);
                
        %%% pixel correlation %%%
        tuset = vld_prd_slct(datasmthf(i, :), cutofff(i), pkcutofff(i));
        tuse = false(1, nf);
        tuse(round(gsw / 2): round(gsw / 2) + length(tuset) - 1) = tuset;

        %%% coarse spatial init %%%
        saa = reshape(sa, lsml, nf);
        daa = saa(:, tuse);
        datt = datusef(i, :);
        dat = datt(tuse);
        ind = corr(daa', dat') > cthres;
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

        if mod(i, round(nseed / 10)) == 0
            disp(['Done #', num2str(i), '/', num2str(nseed)])
        end
    end
    
    %%% get background roi and sig %%%
    %%% compute residual spatial & temporal %%%
    roiid = find(max(roi, [], 2) > 0);
    tominus = zeros(1, nf);
    for i = 1: length(roiid)
        tmp = roi(roiid(i), :);
        idt = tmp > 0;
        tominus = tominus + max(tmp(idt)' .* sig(idt, :), [], 1);
    end
    bgf = (imeantf * pixh * pixw - tominus) / (pixh * pixw);
    bg = (reshape(imeanf * nf, pixh * pixw, 1) - max(roi .* sum(sig, 2)', [], 2)) / nf;

    disp('Done manual pix select')
end