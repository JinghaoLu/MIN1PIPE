function [Y, mxfn, xff, ldf, idf] = inter_section(Y, stt, stp, pixs, scl, sigma_x, sigma_f, sigma_d)
% register frames across stable sections
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
    
    %% compute feature maps %%
    %%% find the feature map of each stable section %%%
    mxall = feature1_comp(Y, stt, stp);
    mxall = normalize(mxall);
    
    %%% 2nd feature map (mxalla) of the 1st (mxall): enhance signal pixels %%%
    mag = 40;
    flag = 0; %%% no anidenoise %%%
    mxalla = feature2_comp(mxall, flag, mag);
    mxalla = normalize(mxalla);
    
    %% LK-LogDemons hierarchical registration %%
    [mxfn, xff, ldf, idf] = lk_ld_hier(mxall, mxalla, pixs, scl, sigma_x, sigma_f, sigma_d);
    
    %% warp all the stable section frames %%
    Yuse = cell(1, length(stt));
    for i = 1: length(stt)
        if i == 1
            Yuse{i} = Y(:, :, 1: stt(i + 1) - 1);
        elseif i == length(stt)
            Yuse{i} = Y(:, :, stt(i): end);
        else
            Yuse{i} = Y(:, :, stt(i): stt(i + 1) - 1);
        end
    end
    Yuse = logdemons_warp_layers(Yuse, xff, ldf);
    for i = 1: length(stt)
        if i == 1
            Y(:, :, 1: stt(i + 1) - 1) = Yuse{i};
        elseif i == length(stt)
            Y(:, :, stt(i): end) = Yuse{i};
        else
            Y(:, :, stt(i): stt(i + 1) - 1) = Yuse{i};
        end
    end
end