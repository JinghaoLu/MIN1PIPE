function [m_out, imaxf, overwrite_flag, imx2, imn2, ibmean] = neural_enhance(m_in, filename, Params)
% batch version of anisotropic diffusion & background removal
%   Jinghao Lu, 07/01/2018

    h = tic;
    msg = 'Overwrite neural enhanced & registrated .mat file (data)? (y/n)';
    overwrite_flag = judge_file(filename, msg);
    
    %% batch preparation %%
    idt = strfind(filename, '_reg.mat');
    fname = filename(1: idt - 1);
    [pixh, pixw, nf] = size(m_in, 'frame_all');
    ttype = class(m_in.frame_all(1, 1, 1));
    stype = parse_type(ttype);
    nsize = pixh * pixw * nf * stype * feature('numCores') + pixh * pixw * 8 * 8 * 2; %%% size of stype * number of parallel cores + 8 convoluton * 8 for double * 2 variables %%%
    nbatch = batch_compute(nsize);
    ebatch = ceil(nf / nbatch);
    idbatch = [1: ebatch: nf, nf + 1];
    nbatch = length(idbatch) - 1;
    
    if overwrite_flag
        if exist(filename, 'file')
            delete(filename)
        end
        
        %% anisotropic diffusion preparation %%
        isparaad = Params.anidenoise_ispara;
        szad = Params.neuron_size;
        iter = Params.anidenoise_iter;
        dt = Params.anidenoise_dt;
        kappa = Params.anidenoise_kappa;
        opt = Params.anidenoise_opt;
        
        %% background removal preparation %%
        isparabr = Params.bg_remove_ispara;
        szbr = Params.neuron_size;
        
        %% batch neural enhancing %%
        imaxf = zeros(pixh, pixw);
        iminf = zeros(pixh, pixw);
        ibmean = zeros(pixh, pixw);
        for i = 1: nbatch
            %%% get the current batch frames %%%
            tmp = m_in.frame_all(1: pixh, 1: pixw, idbatch(i): idbatch(i + 1) - 1);
            
            %%% clean dirts %%%
            disp(['Begin dirts clean #', num2str(i), '/', num2str(nbatch), ' batch'])
            Ydcln = dirt_clean(tmp, szad, isparaad);
            Ydcln = Ydcln + tmp;
            disp(['Done dirts clean #', num2str(i), '/', num2str(nbatch), ' batch'])
%             clear tmp            
            
            %%% anisotropic diffusion %%%
            disp(['Begin anisotropic diffusion #', num2str(i), '/', num2str(nbatch), ' batch'])
            Yblur = anidenoise(Ydcln, szad, isparaad, iter, dt, kappa, opt);
%             Yblur = anidenoise(tmp, szad, isparaad, iter, dt, kappa, opt);
            disp(['Done anisotropic diffusion #', num2str(i), '/', num2str(nbatch), ' batch'])
            clear Ydcln
%             clear tmp
            
            %%% background removal %%%
            disp(['Begin background removal #', num2str(i), '/', num2str(nbatch), ' batch'])
            reg = bg_remove(Yblur, szbr, isparabr); 
            disp(['Done background removal #', num2str(i), '/', num2str(nbatch), ' batch'])
            clear Yblur
            
            %%% get maxs and mins %%%
            imaxf = max(cat(3, max(reg, [], 3), imaxf), [], 3);
            iminf = min(cat(3, min(reg, [], 3), iminf), [], 3);
            bground = tmp - reg;
            
            %%% remove large uniform background %%%
            parfor ii = 1: size(bground, 3)
                bground(:, :, ii) = imtophat(bground(:, :, ii), strel('disk', round(szad * 6)));
            end
            ibmean = ibmean + sum(bground, 3);
            
            %%% save Ydebg %%%
            disp(['Begin saving Ydebg #', num2str(i), '/', num2str(nbatch), ' batch'])
            savef(filename, 2, 'reg')
            disp(['Done saving Ydebg #', num2str(i), '/', num2str(nbatch), ' batch'])
            clear reg
            clear tmp
        end
        toc(h)
        
        %% normalize Ydebg %%
        ibmean = ibmean / nf;
        imx2 = max(imaxf(:));
        imn2 = min(iminf(:));
        imaxf = normalize(imaxf);
        m_out = normalize_batch(filename, 'reg', imx2, imn2, idbatch);
        save([fname, '_supporting.mat'], 'imx2', 'imn2', 'imaxf', 'ibmean', '-append')
    else
        load([fname, '_supporting.mat'], 'imx2', 'imn2', 'imaxf', 'ibmean')
        m_out = matfile(filename, 'writable', true);
    end
    timeh = toc(h);
    disp(['Done neural enhancing, time ', num2str(timeh)])
end