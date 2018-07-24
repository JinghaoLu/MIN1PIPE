function m = nonstable_section(m, stt, stp, pixs, scl, sigma_x, sigma_f, sigma_d)
% register frames within nonstable sections
%   Jinghao Lu, 02/02/2018

    [pixh, pixw, nf] = size(m, 'reg');
    if nargin < 4 || isempty(pixs)
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
    
    if stp(end) ~= nf
        bstt = [bstt; stp(end) + 1];
        bstp = [bstp; nf];
    end
    
    %% further adjust sections for memory fitness %%
    alen = bstp - bstt;
    nff = max(alen);
    ttype = class(m.reg(1, 1, 1));
    stype = parse_type(ttype);
    nsize = pixh * pixw * nff * stype; %%% size of single %%%
    nbatch = batch_compute(nsize);
    ebatch = round(nff / nbatch);
    idx = find(alen > ebatch);
    for i = 1: length(idx)
        tmp = bstt(idx(i)): ebatch: stp(idx(i));
        tmp = tmp(:);
        bstt = [bstt; round(tmp(2: end - 1))];
        bstp = [bstp; round(tmp(2: end - 1)) - 1];
    end
    bstt = sort(bstt);
    bstp = sort(bstp);    

    %% compute batch %%
    df = bstp - bstt + 1;
    dfc = cumsum(df);
    nff = dfc(end);
    nsize = pixh * pixw * nff * 8 * feature('numCores'); %%% size of double %%%
    nbatch = batch_compute(nsize);
    ebatch = ceil(nff / nbatch);

    i = 1;
    idbatch = zeros(1, nbatch);
    while dfc(end) > 0
        idtmp = find(dfc < ebatch, 1, 'last');
        idbatch(i) = idtmp;
        dfc = dfc - dfc(idtmp);
        i = i + 1;
    end
    nbatch = i - 1;
    idbatch = [0, idbatch];

    for i = 1: nbatch
        %% data preparation %%
        %%% get the frames of the nonstable sections %%%
        regpara = cell(1, idbatch(i + 1) - idbatch(i));
        for ii = idbatch(i) + 1: idbatch(i + 1)
            regpara{ii - idbatch(i)} = m.reg(1: pixh, 1: pixw, max(1, bstt(ii) - 1): min(nf, bstp(ii) + 1));
        end
        
        %%% augment edge blocks %%%
        if bstt(1) == 1
            if i == 1
                regpara{1} = cat(3, regpara{1}(:, :, end), regpara{1});
            end
        end
        
        if bstp(end) == nf
            if i == nbatch
                regpara{end} = cat(3, regpara{end}, regpara{end}(:, :, 1));
            end
        end
        
        %% LK-LogDemons %%
        stof = idbatch(i);
        parfor ii = idbatch(i) + 1: idbatch(i + 1)
            regcur = double(regpara{ii - stof});
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
            for j = 1: ncur - 2
                xfcur{j} = [xfcur{j}, xft];
                ldcur{j} = [ldcur{j}, {cellfun(@(x, y) cat(3, x, y), sxt, syt, 'uniformoutput', false)}];
            end
            
            %%% warp the current nonstable section %%%
            regcur = logdemons_warp_layers(squeeze(mat2cell(regcur(:, :, 2: end - 1), pixh, pixw, ones(1, ncur - 2))), xfcur, ldcur);
            regcur = reshape(cell2mat(regcur(:)'), pixh, pixw, ncur - 2);
            regpara{ii - stof} = regcur;
            if length(bstt) < 10
                disp(['Done nonstable-LogDemons section # ', num2str(ii), '/', num2str(length(bstt))])
            else
                if mod(ii, round(length(stt) / 10)) == 0
                    disp(['Done nonstable-LogDemons section # ', num2str(ii), '/', num2str(length(bstt))])
                end
            end
        end
        
        %%% combine nonstable sections into the full tensor %%%
        for ii = idbatch(i) + 1: idbatch(i + 1)
            m.reg(1: pixh, 1: pixw, bstt(ii): bstp(ii)) = regpara{ii - idbatch(i)};
        end
    end
end