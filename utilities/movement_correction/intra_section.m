function m = intra_section(m, stt, stp, pixs, scl, sigma_x, sigma_f, sigma_d, flag, maskc)
% register frames within stable sections
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
    
    if nargin < 9 || isempty(flag)
        flag = 1;
    end

    if nargin < 10 || isempty(maskc)
        maskc = true(pixh, pixw);
    end
    
    %%% if real selection, jump the logdemons part %%%
    if flag
        sclld = scl;
    else
        sclld = scl;
    end
        
    %% prepare for parallel computing %%
    df = stp - stt + 1;
    dfc = cumsum(df);
    nff = dfc(end);
    %     ttype = class(m.reg(1, 1, 1));
    %     stype = parse_type(ttype);
    %     nsize = pixh * pixw * nff * stype; %%% size of single %%%
    nsize = pixh * pixw * nff * 8; %%% size of single %%%
    nbatch = batch_compute(nsize);
    ebatch = ceil(nff / nbatch);
    
    i = 1;
    idbatch = zeros(1, nbatch);
    while dfc(end) > 0
        idtmp = find(dfc <= ebatch, 1, 'last');
        idbatch(i) = idtmp;
        dfc = dfc - dfc(idtmp);
        i = i + 1;
    end
    nbatch = i - 1;
    idbatch = [0, idbatch];
    
    %% intrasection registration %%
    for i = 1: nbatch
        regtpara = cell(1, idbatch(i + 1) - idbatch(i));
        for j = idbatch(i) + 1: idbatch(i + 1)
            regtpara{j - idbatch(i)} = m.reg(1: pixh, 1: pixw, stt(j): stp(j));
        end
        
        %%% find a best registration for frames within each stable section %%%
        stof = idbatch(i);
        parfor ii = 1: length(regtpara)
            regtcur = double(regtpara{ii});
            ncur = size(regtcur, 3);
            [mxcur, xfcur, ldcur, idcur] = lk_ld_hier(regtcur, [], pixs, scl, sigma_x, sigma_f, sigma_d, sclld, maskc); %%% lk_loop+lk_cluster+logdemons_loop %%%
            regtcur = logdemons_warp_layers(squeeze(mat2cell(regtcur, pixh, pixw, ones(1, ncur))), xfcur, ldcur);
            regtcur = reshape(cell2mat(regtcur(:)'), pixh, pixw, ncur);
            regtpara{ii} = regtcur;
            
            %%%% release worker memory %%%%
            regtcur = [];
            ldcur = [];
            
            if length(stt) < 10
                disp(['Done intra_section # ', num2str(ii + stof), '/', num2str(length(stt)), ', batch ', num2str(i), '/', num2str(nbatch)])
            else
                if mod(ii, round(length(stt) / 10)) == 0
                    disp(['Done intra_section # ', num2str(ii + stof), '/', num2str(length(stt)) ', batch ', num2str(i), '/', num2str(nbatch)])
                end
            end
        end
        
        %%% write stable sections to the file %%%
        for j = idbatch(i) + 1: idbatch(i + 1)
            m.reg(1: pixh, 1: pixw, stt(j): stp(j)) = regtpara{j - idbatch(i)};
        end
    end
end