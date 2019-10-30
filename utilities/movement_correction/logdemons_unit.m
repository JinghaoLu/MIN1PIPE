function [imgo, sxfn, syfn] = logdemons_unit(imref, imcur, pixs, scl, sigma_x, sigma_f, sigma_d, maskc)
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

    if nargin < 8 || isempty(maskc)
        maskc = true(size(imref));
    end
    
    %% initialization %%
    pixthres = scl * pixs;
    if gpuDeviceCount >=1
        isgpu = true;
        try
            gpuDevice(1);
        catch
            isgpu = false;
        end
    else
        isgpu = false;
    end
    nlevel = 3;
    maxiter = 20;
    mq = 0.001;
    demr = 255;
    sxfn = {};
    syfn = {};
    sco = get_trans_score(cat(3, imref, imcur), [], [], [], mq, maskc);
    flag = true;
    dsc = 1;
    count = 1;
        
    %% LogDemons registration %%
    while (sco > pixthres && dsc > 0) || flag
        [imcurt, sx, sy] = logdemons(demr * imref, demr * imcur, isgpu, nlevel, sigma_x, sigma_f, sigma_d); %%% with sigma_x = 1 %%%
        imcurt = imcurt / demr;
        sc = get_trans_score(cat(3, imref, gather(imcurt)), [], [], [], mq, maskc);
        dsc = sco - sc;
        flag = false;
        if dsc > 0
            imcur = imcurt;
            if isa(sx, 'gpuArray')
                sxfn{count} = gather(sx);
                syfn{count} = gather(sy);
            else
                sxfn{count} = sx;
                syfn{count} = sy;
            end
            count = count + 1;
            sco = sc;
        else
            if sigma_x > 1
                sigma_x = sigma_x / 2;
                flag = true;
            end
        end
        
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