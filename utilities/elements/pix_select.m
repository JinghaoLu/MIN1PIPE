function [roi, sig, bg, bgf, idusef, datasmthf, cutofff, pkcutofff] = pix_select(m, mask, sz, Fs, sigthres, corrthres)
% select real seeds for downstream processing
%   Jinghao Lu 05/20/2016

    hpix = tic;
    %% initialization %%
    %%% initialize parameters %%%
    if nargin < 3 || isempty(sz)
        defpar = default_parameters;
        sz = defpar.neuron_size;
    end
    
    if nargin < 4 || isempty(Fs)
        defpar = default_parameters;
        Fs = defpar.Fsi_new;
    end
    
    if nargin < 5 || isempty(sigthres)
        defpar = default_parameters;
        sigthres = defpar.pix_select_sigthres; %%% signal salience thres %%%
    end
    
    if nargin < 6 || isempty(corrthres)
        defpar = default_parameters;
        corrthres = defpar.pix_select_corrthres;
    end

    %%% prepare file path %%%
    pathname = mfilename('fullpath');
    mns = mfilename;
    lname = length(mns);
    pathname = pathname(1: end - lname);
    
    %% over-complete set of seeds selection %%
    [pixh, pixw, nf] = size(m, 'reg');
    
    %%% random frames median-max projection %%%
    niter = 10; 
    nselmax = 500; %%% each time only a small portion is collected %%%
    nsel = min(floor(nf / niter), nselmax);
    stype = parse_type(class(m.reg(1, 1, 1)));
    nsize = pixh * pixw * nsel * stype; %%% size of single %%%
    nbatch = batch_compute(nsize);
    ebatch = ceil(nsel / nbatch);
    idbatch = [1: ebatch: nsel, nsel + 1];
    nbatch = length(idbatch) - 1;
    mxmx = zeros(pixh, pixw, niter);
    maxall = zeros(pixh, pixw);
    imeantf = zeros(1, nf);
    imeanf = zeros(pixh, pixw);
    for i = 1: niter
        istt = (i - 1) * nsel;
        maxt = zeros(pixh, pixw);
        for j = 1: nbatch
            tmp = m.reg(1: pixh, 1: pixw, istt + idbatch(j): istt + idbatch(j + 1) - 1);
            maxt = max(max(tmp, [], 3), maxt);
            maxall = max(maxall, maxt);
            imeanf = (imeanf * (istt + idbatch(j) - idbatch(1)) + sum(tmp, 3)) / (istt + idbatch(j + 1) - idbatch(1));
            imeantf(istt + idbatch(j): istt + idbatch(j + 1) - 1) = mean(reshape(tmp, pixh * pixw, idbatch(j + 1) - idbatch(j)), 1);
        end
        mxmx(:, :, i) = imregionalmax(maxt);
%         disp(num2str(i))
    end
    
    nff = nf - nsel * niter;
    if nff > 0
        nsize = pixh * pixw * nff * stype; %%% size of single %%%
        nbatch = batch_compute(nsize);
        ebatch = ceil(nff / nbatch);
        idbatch = [1: ebatch: nff, nff + 1];
        nbatch = length(idbatch) - 1;
        for i = 1: nbatch
            tmp = m.reg(1: pixh, 1: pixw, nsel * niter + idbatch(i): nsel * niter + idbatch(i + 1) - 1);
            maxall = double(max(max(tmp, [], 3), maxall));
            imeanf = (imeanf * (nsel * niter + idbatch(i) - idbatch(1)) + double(sum(tmp, 3))) / (nsel * niter + idbatch(i + 1) - idbatch(1));
            imeantf(nsel * niter + idbatch(i): nsel * niter + idbatch(i + 1) - 1) = double(mean(reshape(tmp, pixh * pixw, idbatch(i + 1) - idbatch(i)), 1));
        end
    end
    
    proj1 = (sum(mxmx, 3) > 0) | (imregionalmax(imgaussfilt(maxall, 1)));
    maskc = bwconvhull(maxall > intensity_filter(maxall));
    if sum(maskc(:)) == 0
        maskc = true(size(maxall));
    end
    disp('Done randomized seeds init')
    toc(hpix)
    
    %% coarse seeds refinement: GMM %%
%     mx = proj1 .* mask .* maskc;
    mx = proj1;
    mxs = find(reshape(mx, pixh * pixw, 1));
    datause = zeros(length(mxs), nf);
    stype = parse_type(class(m.reg(1, 1, 1)));
    nsize = pixh * pixw * nf * stype; %%% size of single %%%
    nbatch = batch_compute(nsize);
    ebatch = ceil(nf / nbatch);
    idbatch = [1: ebatch: nf, nf + 1];
    nbatch = length(idbatch) - 1;
    for i = 1: nbatch
        tmp = m.reg(1: pixh, 1: pixw, idbatch(i): idbatch(i + 1) - 1);
        tmp = reshape(tmp, pixh * pixw, idbatch(i + 1) - idbatch(i));
        datause(:, idbatch(i): idbatch(i + 1) - 1) = tmp(mxs, :);
%         disp(num2str(i))
    end
    
    %%% test direct gmm %%%
    dtmx = prctile(datause', 99.9)';
    dtmn = prctile(datause', 0.1)';
    scmn = dtmx - dtmn;
    G = fitgmdist(scmn, 2);
    [~, iuse] = max(G.mu);
    scall = posterior(G, scmn);
    idxuse2 = scall(:, iuse) > 0.5;
    idxuse3 = scmn > intensity_filter(maxall);
    mxuse = mxs(idxuse2 & idxuse3);
    datause2 = datause(idxuse2 & idxuse3, :);
    time = toc(hpix);
    disp(['Done coarse selection, ', num2str(time), ' seconds'])
    
    %% signal level %%    
    %%% max intensity filter %%%
    ithres = intensity_filter(maxall);
    [~, idm, ~] = intersect(mxuse, find(maxall(:) > ithres));
    datuse = datause2(idm, :);
    iduse = mxuse(idm);
    
    %% neuronal shape filter %%
    krnt = fspecial('gaussian', [pixh, pixw], floor(sz / 2));
    krn = krnt(floor(pixh / 2 - sz): ceil(pixh / 2 + sz), floor(pixw / 2 - sz): ceil(pixw / 2 + sz));
    maxcorr = normxcorr2(krn, maxall);
    [hk, wk] = size(krn);
    maxcorr = maxcorr(ceil(hk / 2): end - ceil(hk / 2) + 1, ceil(wk / 2): end - ceil(wk / 2) + 1);
    [~, idm, ~] = intersect(iduse, find(maxcorr(:) > 0));
    datuse = datuse(idm, :);
    iduse = iduse(idm);    
        
    %% seeds refine: preparation for temporal calcium spike dynamics %%
    gsw = 1 * Fs;
    datasmth = convn(datuse, gausswin(gsw)' / sum(gausswin(gsw)), 'valid');
    st = sort(datasmth, 2);
    bs = zeros(size(st));
    for i = 1: size(st, 1)
        bs(i, :) = linspace(st(i, 1), st(i, end), size(st, 2));
    end
    stcr = st - bs;
    [~, minp] = min(stcr, [], 2);
    [~, maxp] = max(stcr, [], 2);
    minids = sub2ind(size(st), (1:size(st,1))', minp);
    maxids = sub2ind(size(st), (1:size(st,1))', maxp);
    cutoff = (st(minids) + st(maxids)) / 2;
    pkcutoff = st(minids);

    %% seeds refinement: exclude seeds with no calcium component %%
    ksscr = zeros(1, length(iduse));
    parfor i = 1: length(iduse)
        tmp = datuse(i, :);
        ksscr(i) = kstest(zscore(tmp));
    end
    ksscr = logical(ksscr);
    datuse = datuse(ksscr, :);
    datasmth = datasmth(ksscr, :);
    cutoff = cutoff(ksscr);
    pkcutoff = pkcutoff(ksscr);
    iduse = iduse(ksscr);

    %% seeds refinement: exclude non calcium property %%
    %%% skewness %%%
    skw = skewness(datuse');
    idd = find(skw > 0);
    datuse = datuse(idd, :);
    datasmth = datasmth(idd, :);
    cutoff = cutoff(idd);
    pkcutoff = pkcutoff(idd);
    iduse = iduse(idd);            
    time = toc(hpix);
    disp(['Done refinement, ', num2str(time), ' seconds'])    
    
    %% seeds correlating %%   
    [idusef, datasmthf, datusef, cutofff, pkcutofff] = seeds_merge(maxall, iduse, datuse, datasmth, cutoff, pkcutoff, sz, corrthres);
    
    %% seeds clean %%
    maskall = mask .* maskc;
    idinmask = maskall(idusef) > 0;
    idusef = idusef(idinmask);
    datasmthf = datasmthf(idinmask, :);
    datusef = datusef(idinmask, :);
    cutofff = cutofff(idinmask);
    pkcutofff = pkcutofff(idinmask);
        
    %% rnn classifier %%
%     idend = find(pathname == filesep, 2, 'last');
%     load([pathname(1: idend(1)), 'RNN', filesep, 'rnn_lstm_model.mat'])
%     [lbl, ~] = seeds_cleansing_rnn(datuse, model);
%     lbl = find(lbl == 1);
%     datuse = datuse(lbl, :);
%     datasmth = datasmth(lbl, :);
%     cutoff = cutoff(lbl);
%     pkcutoff = pkcutoff(lbl);
%     iduse = iduse(lbl);
%     disp('Done LSTM')
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     %% seeds correlating: post %%   
%     [idusef, datasmthf, datusef, cutofff, pkcutofff] = seeds_merge(frame, sz, iduset, datuset, datasmtht, cutofft, pkcutofft);

    %% spatiotemporal initialization %%
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
        sa = m.reg(rg(1, 1): rg(1, 2), rg(2, 1): rg(2, 2), 1: nf);
                
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
    
    time = toc(hpix);
    disp(['Done pix select, time: ', num2str(time), ' seconds'])
end

