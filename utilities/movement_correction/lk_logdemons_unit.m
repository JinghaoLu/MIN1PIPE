function [tmp, xff, sxx, syy] = lk_logdemons_unit(tmp, se, pixs, scl, sigma_x, sigma_f, sigma_d, flag, maskc) 
% LK logdemons registration unit
%   Jinghao Lu, 10/20/2018

    %%% initialize %%%
    [pixh, pixw, nn] = size(tmp);
    
    if nargin < 2 || isempty(se)
        defpar = default_parameters;
        se = defpar.neuron_size;
    end
    
    if nargin < 3 || isempty(pixs)
        pixs = min(pixh, pixw);
    end
    
    if nargin < 4 || isempty(scl)
        defpar = default_parameters;
        scl = defpar.mc_scl;
    end
    
    if nargin < 5 || isempty(sigma_x)
        defpar = default_parameters;
        sigma_x = defpar.mc_sigma_x;
    end
    
    if nargin < 6 || isempty(sigma_f)
        defpar = default_parameters;
        sigma_f = defpar.mc_sigma_f;
    end
    
    if nargin < 7 || isempty(sigma_d)
        defpar = default_parameters;
        sigma_d = defpar.mc_sigma_d;
    end
    
    if nargin < 8 || isempty(flag)
        flag = false;
    end
    
    if nargin < 9 || isempty(maskc)
        maskc = true(pixh, pixw);
    end
    
    xff = cell(1, nn);
    sxx = cell(1, nn);
    syy = cell(1, nn);
    p0 = [0; 0];

    %%% imgaussfilt %%%
    for j = 1: nn
        tmp(:, :, j) = imgaussfilt(tmp(:, :, j), 1);
    end
    
    %%% first LK track %%%
    PNs = zeros(2, nn - 1);
    for j = 1: nn - 1
        imref = tmp(:, :, j);
        imcur = tmp(:, :, j + 1);
        [PNt, imgt, sco] = lk_ref_track(imcur, imref, maskc);
        scc = get_trans_score_ref(imgt, imref, maskc);
        if scc < sco
            PN = PNt;
        else
            PN = p0;
        end
        PNs(:, j) = PN;
    end
    PNs = cumsum(PNs, 2);
    
    %%% warp LK %%%
    xff{1} = p0;
    for j = 1: nn - 1
        tmp(:, :, j + 1) = lk2_warp(tmp(:, :, j + 1), PNs(:, j), p0);
        xff{j + 1} = PNs(:, j);
    end
    
    %%% feature map %%%
    tmpp = tmp;
    if flag
        for j = 1: nn
            tmpp(:, :, j) = feature2_comp(tmpp(:, :, j), [], 100);
        end
    end
    
    %%% 1st logdemons track %%%
    sxx{1} = {zeros(pixh, pixw)};
    syy{1} = {zeros(pixh, pixw)};
    for j = 1: nn - 1
        imref = max(tmpp(:, :, 1: j), [], 3);
        imcur = tmpp(:, :, j + 1);
        
        [img, sxt, syt] = logdemons_unit(imref, imcur, pixs, scl, sigma_x, sigma_f, sigma_f, maskc);
        if isempty(sxt)
            sxt = {zeros(pixh, pixw)};
            syt = {zeros(pixh, pixw)};
        end
        
        tmpp(:, :, j + 1) = img;
        for k = 1: length(sxt)
            tmp(:, :, j + 1) = iminterpolate(tmp(:, :, j + 1), sxt{k}, syt{k});
        end
        
        sxx{j + 1} = sxt;
        syy{j + 1} = syt;
    end
    
    %%% 1st imopen %%%
    for j = 1: nn
        tmp(:, :, j) = tmp(:, :, j) - imopen(tmp(:, :, j), strel('disk', se));
    end
    
    %%% 2nd logdemons track %%%
    imref = max(tmp, [], 3);
    imref = imref - imopen(imref, strel('disk', se));
    for j = 1: nn
        imcur = tmp(:, :, j);
        
        [img, sxt, syt] = logdemons_unit(imref, imcur, pixs, scl, sigma_x, sigma_f, sigma_d, maskc);
        if isempty(sxt)
            sxt = {zeros(pixh, pixw)};
            syt = {zeros(pixh, pixw)};
        end
        
        tmp(:, :, j) = img;
        
        sxx{j} = [sxx{j}, sxt];
        syy{j} = [syy{j}, syt];
    end
    
    %%% 2nd imopen %%%
    for j = 1: nn
        tmp(:, :, j) = tmp(:, :, j) - imopen(tmp(:, :, j), strel('disk', se));
    end
end