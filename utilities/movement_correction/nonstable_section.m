function Y = nonstable_section(Y, stt, stp, pixs, scl, sigma_x, sigma_f, sigma_d)
% register frames within nonstable sections
%   Jinghao Lu, 02/02/2018

    if nargin < 4 || isempty(pixs)
        [pixh, pixw, ~] = size(Y);
        pixs = min(pixh, pixw);
    end
    
    if nargin < 5 || isempty(scl)
        defpar = default_parameters;
        scl = defpar.mc_scl;
    end
    
    if nargin < 6 || isempty(sigma_x)
        defpar = default_parameters;
        sigma_x = defpar.mc_sigma_x;
    end
    
    if nargin < 7 || isempty(sigma_f)
        defpar = default_parameters;
        sigma_f = defpar.mc_sigma_f;
    end
    
    if nargin < 8 || isempty(sigma_d)
        defpar = default_parameters;
        sigma_d = defpar.mc_sigma_d;
    end
    
    %% generate nonstable clusters %%
    %%% compute the range of the nonstable sections %%%
    [pixh, pixw, nframes] = size(Y);
    bstt = stp(1: end - 1) + 1;
    bstp = stt(2: end) - 1;
    
    tmp = bstp - bstt;
    id = tmp >= 0;
    bstt = bstt(id);
    bstp = bstp(id);
    
    if stt(1) ~= 1
        bstt = [1; bstt];
        bstp = [stt(1) - 1; bstp];
    end
    
    if stp(end) ~= nframes
        bstt = [bstt; stp(end) + 1];
        bstp = [bstp; nframes];
    end
    
    %%% get the frames of the nonstable sections %%%
    regpara = cell(1, length(bstt));
    for i = 1: length(bstt)
        regpara{i} = Y(:, :, max(1, bstt(i) - 1): min(nframes, bstp(i) + 1));
    end
    
    %%% augment edge blocks %%%
    if bstt(1) == 1
        regpara{1} = cat(3, regpara{1}(:, :, end), regpara{1});
    end
    if bstp(end) == nframes
        regpara{end} = cat(3, regpara{end}, regpara{end}(:, :, 1));
    end
    
    %% LK-LogDemons %%
    parfor ii = 1: length(bstt)
        regcur = regpara{ii};
        regcurt = normalize(regcur);
        ncur = size(regcur, 3);
                
        %%% register with LK-LogDemons hierarchical framework %%%
        [mxcur, xfcur, ldcur, idcur] = lk_ld_hier(regcurt(:, :, 2: end - 1), [], pixs, scl, sigma_x, sigma_f, sigma_d); %%% simply register pure nonstable frames, then corrected to the previous stable section %%%
        
        %%% register back to the first frame / previous stable section %%%
        imref = regcurt(:, :, 1);
        imcur = mxcur(:, :, idcur{1} == 1);
        [scrtc, imgt, xft] = klt_ref_track(imcur, imref, 20); %%% run more times of KLT tracking to ensure best result %%%
        scrto = get_trans_score_ref(imcur, imref);
        if scrtc < scrto
            imcur = imgt;
        else
            xft = {affine2d(diag(ones(1, 3)))};
        end
        [imgo, sxt, syt] = logdemons_unit(imref, imcur, pixs, scl, sigma_x, sigma_f, sigma_d);
        
        %%% convert final LK-LogDemons transforms %%%
        for i = 1: ncur - 2
            xfcur{i} = [xfcur{i}, xft];
            ldcur{i} = [ldcur{i}, {cellfun(@(x, y) cat(3, x, y), sxt, syt, 'uniformoutput', false)}];
        end
        
        %%% warp the current nonstable section %%%
        regcur = logdemons_warp_layers(squeeze(mat2cell(regcur(:, :, 2: end - 1), pixh, pixw, ones(1, ncur - 2))), xfcur, ldcur);
        regcur = reshape(cell2mat(regcur(:)'), pixh, pixw, ncur - 2);
        regpara{ii} = regcur;
        if length(bstt) < 10
            disp(['Done nonstable-LogDemons section # ', num2str(ii), '/', num2str(length(bstt))])
        else            
            if mod(ii, round(length(stt) / 10)) == 0
                disp(['Done nonstable-LogDemons section # ', num2str(ii), '/', num2str(length(bstt))])
            end
        end
    end
    
    %%% combine nonstable sections into the full tensor %%%
    for i = 1: length(bstt)
        Y(:, :, bstt(i): bstp(i)) = regpara{i};
    end
end