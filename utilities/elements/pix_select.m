function [roi, sig, bg, bgf, idusef, datasmthf, cutofff, pkcutofff] = pix_select(frame, mask, sz, Fs, sigthres, corrthres)
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
    [pixh, pixw, nf] = size(frame);
    
    %%% random frames median-max projection %%%
    niter = 10; 
    nselmax = 500; %%% each time only a small portion is collected %%%
    nsel = min(floor(nf / niter), nselmax);
    nfuse = niter * nsel;
    datap = frame(:, :, 1: nfuse);
    datap = reshape(datap, pixh, pixw, niter, nsel);
    datap = max(datap, [], 4);
    mxmx = zeros(pixh, pixw, niter);
    for i = 1: niter
        datamxt = datap(:, :, i);
        mxmx(:, :, i) = imregionalmax(datamxt);
%         disp(num2str(i))
    end
    proj1 = (sum(mxmx, 3) > 0) | (imregionalmax(imgaussfilt(max(frame, [], 3), 1)));
    disp('Done randomized seeds init')
    toc(hpix)
    
    %% coarse seeds refinement: GMM %%
    mx = proj1 .* mask;
    mxs = find(reshape(mx, pixh * pixw, 1));
    datause = reshape(frame, pixh * pixw, nf);
    datause = datause(mxs, :);
    
    %%% test direct gmm %%%
    dtmx = prctile(datause', 99.9)';
    dtmn = prctile(datause', 0.1)';
    scmn = dtmx - dtmn;
    G = fitgmdist(scmn, 2);
    [~, iuse] = max(G.mu);
    scall = posterior(G, scmn);
    idxuse2 = scall(:, iuse) > 0.5;
    mxuse = mxs(idxuse2);
    datause2 = datause(idxuse2, :);
    time = toc(hpix);
    disp(['Done coarse selection, ', num2str(time), ' seconds'])
    
    %% signal level %%
    %%% signal level fiter (compared to noise) %%%
    ffts = fft(datause2, [], 2);
    nft = size(ffts, 2);
    nf4 = round(nft / 4);
    rgs = [floor(nft / 2 + 1) - nf4, ceil(nft / 2 + 1) + nf4];
    ffts(:, [1: rgs(1), rgs(2): end]) = 0;
    datt = ifft(ffts, [], 2);
    rgnoise = max(datt(:, 2: end - 1), [], 2) - min(datt(:, 2: end - 1), [], 2);
    rgsig = max(datause2, [], 2) - min(datause2, [], 2);
    idm = rgsig ./ rgnoise > sigthres;
    datause2 = datause2(idm, :);
    mxuse = mxuse(idm);
    
    %%% max intensity filter %%%
    imgmax = max(frame, [], 3);
    nbins = round(pixh * pixw / 10);
    [tmp1, ctrs] = hist(imgmax(:), nbins);
    [~, idm] = max(tmp1);
    idthres = 2 * idm;
    ithres = ctrs(idthres);
    [~, idm, ~] = intersect(mxuse, find(imgmax(:) > ithres));
    datuse = datause2(idm, :);
    iduse = mxuse(idm);
        
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
    time = toc(hpix);
    disp(['Done refinement, ', num2str(time), ' seconds'])

    %% seeds correlating %%   
    [idusef, datasmthf, datusef, cutofff, pkcutofff] = seeds_merge(frame, iduse, datuse, datasmth, cutoff, pkcutoff, sz, corrthres);
        
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
    res = frame;
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
        a = res(rg(1, 1): rg(1, 2), rg(2, 1): rg(2, 2), :);
        sa = frame(rg(1, 1): rg(1, 2), rg(2, 1): rg(2, 2), :);
                
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

    time = toc(hpix);
    disp(['Done pix select, time: ', num2str(time), ' seconds'])
end

