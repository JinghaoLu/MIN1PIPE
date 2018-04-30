function [imgo, sxfn, syfn] = logdemons_unit(imref, imcur, pixs, scl, sigma_x, sigma_f, sigma_d)
% LogDemons algorithm in loop
%   Jinghao Lu, 02/22/2018

    if nargin < 3 || isempty(pixs)
        [pixh, pixw] = size(imref);
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

    %% initialization %%
    demthres = 0.01;
    pixthres = scl * pixs;
    isgpu = true;
    nlevel = 3;
    maxiter = 10;
    demr = 255;
    sxfn = {};
    syfn = {};
    sco = get_trans_score(cat(3, imref, imcur));
    dsc = 1;
    flag = true;
    count = 1;
        
    %% LogDemons registration %%
    while (sco > pixthres && dsc > demthres) || flag
        [imcur, sx, sy] = logdemons(demr * imref, demr * imcur, isgpu, nlevel, sigma_x, sigma_f, sigma_d); %%% with sigma_x = 1 %%%
        imcur = imcur / demr;
        if isgpu
            sxfn{count} = gather(sx);
            syfn{count} = gather(sy);
        else
            sxfn{count} = sx;
            syfn{count} = sy;
        end
        sc = get_trans_score(cat(3, imref, gather(imcur)));
        dsc = sco - sc;
        sco = sc;
        count = count + 1;
        flag = false;
        
        if count > maxiter
            break
        end
    end
    
    if isgpu
        imgo = gather(imcur);
    else
        imgo = imcur;
    end
end