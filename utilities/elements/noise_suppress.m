function [m_out] = noise_suppress(m_in, maxall)
% batch version of anisotropic diffusion & background removal
%   Jinghao Lu, 11/18/2019

    h = tic;
    
    %% first intensity filter %%
    [pixh, pixw, nf] = size(m_in, 'reg');
    ithres = intensity_filter(maxall);
    mask = maxall > ithres;
    imax = mask;
    imin = mask;

    %% batch preparation %%
    ttype = class(m_in.reg(1, 1, 1));
    stype = parse_type(ttype);
    nsize = pixh * pixw * nf * stype * feature('numCores'); 
    nbatch = batch_compute(nsize);
    ebatch = ceil(pixw / nbatch);
    idbatch = [1: ebatch: pixw, pixw + 1];
    nbatch = length(idbatch) - 1;
    
    for ib = 1: nbatch
        data = m_in.reg(:, idbatch(ib): idbatch(ib + 1) - 1, :);
        data = reshape(data, [], nf);
        iduse = mask(:, idbatch(ib): idbatch(ib + 1) - 1);
        iduse = find(iduse);
        datause = data(iduse, :);
        
        %%% second kstest %%%
        ksr = zeros(length(iduse), 1);
        parfor i = 1: length(iduse)
            tmp = datause(i, :);
            ksr(i) = kstest(zscore(tmp));
        end
        idm1 = find(ksr > 0);
        data1 = datause(idm1, :);
        
        %%% third gaussian fit %%%
        thress = zeros(length(idm1), 1);
        dthres = 0.01;
        id1 = true(length(idm1), 1);
%         mn = zeros(length(idm1), 1);
        parfor i = 1: length(idm1)
            tmp = data1(i, :);
            [n, edges] = histcounts(tmp);
            bin = (edges(1: end - 1) + edges(2: end)) / 2;
            g = fit(double(bin)', double(n)', 'gauss1', 'lower', [0, -1, 0], 'upper', [1e8, 1, 1]);
%             mn(i) = g.b1;
            y = feval(g, bin);
            xx = n / sum(n);
            yy = y / sum(y);
            d = xx(:) - yy(:);
            threst = bin(max(1, find(d(:) > dthres & bin(:) > g.b1, 1) - 1));
            if isempty(threst)
                id1(i) = false;
            else
                thress(i) = threst;
            end
%             if mod(i, 100) == 0
%                 disp(num2str(i))
%             end
        end
        mn = mean(data1, 2);
        % mn = mean(data1, 2) - 2 * mad(data1, 0, 2);
        mn_use = mn(id1);
        thress_use = thress(id1);
        data2 = data1(id1, :);
        
        %%% fourth rectified below-threshold suppress of adaptive noise %%%
        thressnew = thress_use - mn_use;
        for i = 1: length(thress_use)
            tmp = data2(i, :);
            tmp = tmp - mn_use(i);
            tmpt = sigmf(tmp, [100, thressnew(i)]);
            data2(i, :) = tmpt .* tmp;
            data2(i, :) = data2(i, :) - min(data2(i, :));
        end
        
        %%% minus mean/min and rectify %%%
        % datat = max(data - mean(data, 2), 0);
        data = max(data - mean(data, 2), 0);
        data = data .^ 4;
        data(iduse(idm1(id1)), :) = data2;
        data = reshape(data, pixh, [], nf);
        
        %%% update data %%%
        m_in.reg(:, idbatch(ib): idbatch(ib + 1) - 1, :) = data;
        imax = max(data, [], 3);
        imin = min(data, [], 3);
    end
    
    imax = max(imax(:));
    imin = min(imin(:));
    m_out = normalize_batch(m_in.Properties.Source, 'reg', imax, imin, idbatch);
    
    %%% smooth %%%
    ebatch = ceil(nf / nbatch);
    idbatch = [1: ebatch: nf, nf + 1];
    nbatch = length(idbatch) - 1;
    for ib = 1: nbatch
        data = m_in.reg(:, :, idbatch(ib): idbatch(ib + 1) - 1);
        for i = 1: size(data, 3)
            data(:, :, i) = imgaussfilt(data(:, :, i), 1.5);
        end
        m_in.reg(:, :, idbatch(ib): idbatch(ib + 1) - 1) = data;
    end
    
    timeh = toc(h);
    disp(['Done noise suppression, time ', num2str(timeh)])
end