function [m, mxfn, xff, ldf, idf] = inter_section(m, stt, stp, pixs, scl, sigma_x, sigma_f, sigma_d)
% register frames across stable sections
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
    
    %% compute feature maps %%
    %%% compute batch %%%
    df = stp - stt + 1;
    dfc = cumsum(df);
    nff = dfc(end);
    ttype = class(m.reg(1, 1, 1));
    stype = parse_type(ttype);
    nsize = pixh * pixw * nff * stype; %%% size of single %%%
    nbatch = batch_compute(nsize);
    ebatch = round(nff / nbatch);
    
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

    %%% find the feature map of each stable section %%%
    mxall = feature1_comp(m, stt, stp, idbatch);
    mxall = normalize(mxall);
    
    %%% 2nd feature map (mxalla) of the 1st (mxall): enhance signal pixels %%%
    mag = 40;
    flag = 0; %%% no anidenoise %%%
    mxalla = feature2_comp(mxall, flag, mag);
    mxalla = normalize(mxalla);
    
    %% LK-LogDemons hierarchical registration %%
    [mxfn, xff, ldf, idf] = lk_ld_hier(mxall, mxalla, pixs, scl, sigma_x, sigma_f, sigma_d);
    
    %% warp all the stable section frames %%
    for i = 1: nbatch
        Yuse = cell(1, idbatch(i + 1) - idbatch(i));
        for ii = idbatch(i) + 1: idbatch(i + 1)
            if ii == 1
                Yuse{ii - idbatch(i)} = m.reg(1: pixh, 1: pixw, 1: stt(ii + 1) - 1);
            elseif ii == length(stt)
                Yuse{ii - idbatch(i)} = m.reg(1: pixh, 1: pixw, stt(ii): end);
            else
                Yuse{ii - idbatch(i)} = m.reg(1: pixh, 1: pixw, stt(ii): stt(ii + 1) - 1);
            end
        end
        
        Yuse = logdemons_warp_layers(Yuse, xff(idbatch(i) + 1: idbatch(i + 1)), ldf(idbatch(i) + 1: idbatch(i + 1)));
        
        for ii = idbatch(i) + 1: idbatch(i + 1)
            if ii == 1
                m.reg(1: pixh, 1: pixw, 1: stt(ii + 1) - 1) = Yuse{ii - idbatch(i)};
            elseif ii == length(stt)
                m.reg(1: pixh, 1: pixw, stt(ii): end) = Yuse{ii - idbatch(i)};
            else
                m.reg(1: pixh, 1: pixw, stt(ii): stt(ii + 1) - 1) = Yuse{ii - idbatch(i)};
            end
        end
    end
end