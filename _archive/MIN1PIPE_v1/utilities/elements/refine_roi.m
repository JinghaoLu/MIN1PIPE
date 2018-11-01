function [A, b, C, iduse, datasmthf, cutofff, pkcutofff] = refine_roi(m, C, f, Aold, iduse, noise, datasmthf, cutofff, pkcutofff, ispara)
% [A, b, C, iduse] = refine_roi refine roi by basis pursuit denoising
%   modified from E Pnevmatikakis
%   Jinghao Lu 06/10/2016

    hroi = tic;
    if isempty(gcp('nocreate'))
        parpool(feature('numCores'));
    end
    
    %% initialization %%
    %%% initialize parameters %%%
    if nargin < 10
        defpar = default_parameters;
        ispara = defpar.refine_roi_ispara;
    end
    
    [d1, d2, d3] = size(m, 'reg');
    d = d1 * d2;
    nroi = size(Aold, 2);
    nseed = length(iduse);
    swin = 6;

    %%% imgaussfilt for search area construction %%%
    sarea = zeros(size(Aold));
    cthres = 0.1;
    for i = 1: nroi
        sarea(:, i) = reshape(normalize(imgaussfilt(full(reshape(Aold(:, i), d1, d2)), swin)) > cthres, d, 1);
    end
    time = toc(hroi);
    disp(['Done init, use ', num2str(time), ' seconds'])
    
    %% BPDN main %%
    A = zeros(d, nroi);
    ind2use = find(sum(sarea, 2));
%     [h, w] = ind2sub([d1, d2], ind2use);
    npx = length(ind2use);
    ind = cell(1, npx);
    Cuse = cell(1, npx);
    for i = 1: npx
        ind{i} = find(sarea(ind2use(i), :));
        Cuse{i} = [C(ind{i}, :); f];
    end
    
    nsize = d * d3 * 8; %%% size of double %%%
    nbatch = batch_compute(nsize);
    ebatch = ceil(d / nbatch);
    eb = floor(ebatch / d1);
    idbatch = [1: eb: d2, d2 + 1];
    nbatch = length(idbatch) - 1;

    %%% prepare for parallel %%%
    At = A(ind2use, :);
    noiset = noise(ind2use);
    indcount = 0;
    
    if ispara
        for j = 1: nbatch
            datat = m.reg(1: d1, idbatch(j): idbatch(j + 1) - 1, 1: d3);
            datat = reshape(datat, d1 * (idbatch(j + 1) - idbatch(j)), d3);
            idt = ind2use > d1 * (idbatch(j) - 1) & ind2use <= d1 * (idbatch(j + 1) - 1);
            idt = ind2use(idt) - d1 * (idbatch(j) - 1);
            datat = double(datat(idt, :));
            Ctmp = Cuse(indcount + 1: indcount + length(idt));
            Atmp = At(indcount + 1: indcount + length(idt), :);
            noisetmp = noiset(indcount + 1: indcount + length(idt));
            indtmp = ind(indcount + 1: indcount + length(idt));
            parfor i = 1: length(idt)  % estimate spatial components
                [~, ~, a, ~, ~] = lars_regression_noise(datat(i, :)', Ctmp{i}', 1, noisetmp(i) ^ 2 * d3);
                tmp = Atmp(i, :);
                tmp(indtmp{i}) = a(1: max(1, end - 1))';
                Atmp(i, :) = tmp;
%                 disp(num2str(i))
            end
            Cuse(indcount + 1: indcount + length(idt)) = Ctmp;
            At(indcount + 1: indcount + length(idt), :) = Atmp;
            noiset(indcount + 1: indcount + length(idt)) = noisetmp;
            ind(indcount + 1: indcount + length(idt)) = indtmp;
            indcount = indcount + length(idt);
        end
    else
        for i = 1: npx   % estimate spatial components
            [~, ~, a, ~, ~] = lars_regression_noise(squeeze(double(m.reg(h(i), w(i), 1: d3))), Cuse{i}', 1, noiset(i) ^ 2 * d3);
            tmp = At(i, :);
            tmp(ind{i}) = a(1: max(1, end - 1))';
            At(i, :) = tmp;
%             disp(num2str(i))
        end
    end
    
    %%% update full variable %%%
    A(ind2use, :) = At;
    clear At indt Cuset noiset Yt

    if iscell(A)
        A = cell2mat(A);
    end
    A(isnan(A)) = 0;
    A = A(:, 1: nseed); %%% only get required rois (based on seeds) %%%
    time = toc(hroi);
    disp(['Done BPDN, use ', num2str(time), ' seconds'])
    
    %% Postprocessing %%
    %%% medfilt2 %%%
    mfwin = [3, 3];
    Atmp = reshape(full(A), d1, d2, nseed);
    for i = 1: nseed
        Atmp(:, :, i) = medfilt2(Atmp(:, :, i), mfwin);
    end
    Atmp = reshape(Atmp, d, nseed);
    
    %%% energy select %%%
    ecut = 0.2;
    cumeng = max(Atmp, [], 1);
    icut = cumeng * ecut;
    Atmp2 = bsxfun(@gt, Atmp, icut);
    
    %%% biggest bulk %%%
    %%% have to include the seed %%%
    Atmp3 = zeros(size(Atmp2));
    for i = 1: nseed %%% only compute the required roi %%%
        tmp = reshape(Atmp2(:, i), d1, d2);
        [l, ~] = bwlabeln(full(tmp));
        [x, y] = ind2sub([d1, d2], iduse(i));
        ii = l(x, y);
        Atmp3(:, i) = reshape(l == ii, d, 1) .* Atmp(:, i) .* Atmp2(:, i);
    end
        
    %%% imclose %%%
    se = strel('disk', 3);
    parfor i = 1: nseed
        tmp = reshape(Atmp3(:, i), d1, d2);
        tmp = imclose(tmp, se);
        Atmp3(:, i) = reshape(tmp, d, 1);
    end
    A = Atmp3;
    
    %%% delete empty cell %%%
    inddel = find(sum(A) == 0);
    A = A(:, setdiff(1: nseed, inddel));
    C = C(setdiff(1: nseed, inddel), :);
    datasmthf = datasmthf(setdiff(1: nseed, inddel), :);
    cutofff = cutofff(setdiff(1: nseed, inddel));
    pkcutofff = pkcutofff(setdiff(1: nseed, inddel));
    iduse = iduse(setdiff(1: nseed, inddel));
            
    %%% final smooth of roi %%%
    nseed = length(iduse);
    for i = 1: nseed
        tmp = reshape(A(:, i), d1, d2);
        A(:, i) = reshape(imgaussfilt(tmp, 2), d, 1);
    end
    
    %%% retain only main hill %%%
    for i = 1: length(iduse)
        [y, x] = ind2sub([d1, d2], iduse(i));
        dt = reshape(A(:, i), d1, d2);
        rmax = imregionalmax(dt);
        [ymax, xmax] = ind2sub([d1, d2], find(rmax));
        idx = xmax == x;
        idy = ymax == y;
        ymax = ymax(~(idx & idy));
        xmax = xmax(~(idx & idy));
        np = sum(rmax(:));
        if np > 1
            %%% get all np - 1 lines %%%
            vmint = zeros(1, np - 1);
            for j = 1: np - 1
                ltmp = improfile(dt, [x; xmax(j)], [y; ymax(j)]);
                vmint(j) = min(ltmp);
            end
            vmin = min(vmint);
            
            %%% get main hill %%%
            dtb = dt > vmin;
            [l, n] = bwlabeln(dtb, 4);
            id = l(y, x);
            dtuse = dt .* (l == id);
            dtuse = imgaussfilt(dtuse, 2);
            A(:, i) = dtuse(:);
        end
    end
    
    time = toc(hroi);
    disp(['Done postprocessing, use ', num2str(time), ' seconds'])

    %% background update %%
    nsize = d * d3 * 8; %%% size of double %%%
    nbatch = batch_compute(nsize);
    ebatch = ceil(d / nbatch);
    eb = floor(ebatch / d1);
    idbatch = [1: eb: d2, d2 + 1];
    nbatch = length(idbatch) - 1;
    
    b = zeros(d, 1);
    for i = 1: nbatch
        tmp = m.reg(1: d1, idbatch(i): idbatch(i + 1) - 1, 1: d3);
        tmp = double(reshape(tmp, d1 * (idbatch(i + 1) - idbatch(i)), d3));
        Yf = tmp * f';
        b(d1 * (idbatch(i) - 1) + 1: d1 * (idbatch(i + 1) - 1)) = max((Yf - A(d1 * (idbatch(i) - 1) + 1: d1 * (idbatch(i + 1) - 1), :) * (C * f')) / (f * f'), 0);
    end
    time = toc(hroi);
    disp(['Done refine roi, total time: ', num2str(time), ' seconds'])
end
