function regt = intra_section(Y, stt, stp, pixs, scl, sigma_x, sigma_f, sigma_d, flag)
% register frames within stable sections
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
    
    if nargin < 9 || isempty(flag)
        flag = 1;
    end

    %%% if real selection, jump the logdemons part %%%
    if flag
        sclld = 1;
    end
        
    %% prepare for parallel computing %%
    regt = Y;
    [pixh, pixw, ~] = size(Y);
    regtpara = cell(1, length(stt));
    for i = 1: length(stt)
        regtpara{i} = regt(:, :, stt(i): stp(i));
    end
    
    %% find a best registration for frames within each stable section %%
    parfor ii = 1: length(stt)
        regtcur = regtpara{ii};
        regtcura = regtcur;
        ncur = size(regtcura, 3);
        [mxcur, xfcur, ldcur, idcur] = lk_ld_hier(regtcura, [], pixs, scl, sigma_x, sigma_f, sigma_d, sclld); %%% lk_loop+lk_cluster+logdemons_loop %%%
        regtcura = logdemons_warp_layers(squeeze(mat2cell(regtcura, pixh, pixw, ones(1, ncur))), xfcur, ldcur);        
        regtcura = reshape(cell2mat(regtcura(:)'), pixh, pixw, ncur);
        regtpara{ii} = regtcura;
        if length(stt) < 10
            disp(['Done intra_section # ', num2str(ii), '/', num2str(length(stt))])
        else            
            if mod(ii, round(length(stt) / 10)) == 0
                disp(['Done intra_section # ', num2str(ii), '/', num2str(length(stt))])
            end
        end
    end
    for i = 1: length(stt)
        regt(:, :, stt(i): stp(i)) = regtpara{i};
    end
end