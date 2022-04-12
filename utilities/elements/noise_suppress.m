function [m_out] = noise_suppress(m_in, maxall, Fs, nflag, filename)
% batch version of anisotropic diffusion & background removal
%   Jinghao Lu, 11/18/2019

    h = tic;
    
    if nargin > 4 && ~isempty(filename)
        msg = 'Overwrite post registered .mat file (data)? (y/n)';
        overwrite_flag = judge_file(filename, msg);
        
        if overwrite_flag
            if exist(filename, 'file')
                delete(filename)
            end

            idt1 = strfind(filename, '_reg_post.mat');
            fname_old = [filename(1: idt1 - 1), '_reg.mat'];
            status = copyfile(fname_old, filename, 'f');
        end
        m_in = matfile(filename, 'writable', true);
    else
        overwrite_flag = false;
        filename = [];
    end

    if overwrite_flag || nflag == 1
        %% first intensity filter %%
        [pixh, pixw, nf] = size(m_in, 'reg');
        ithres = intensity_filter(maxall);
        mask = maxall > ithres;
        imax = 1;
        imin = 0;
        
        %% batch preparation %%
        ttype = class(m_in.reg(1, 1, 1));
        stype = parse_type(ttype);
        nsize = pixh * pixw * nf * stype;
        nbatch = batch_compute(nsize);
        ebatch = ceil(pixw / nbatch);
        idbatch = [1: ebatch: pixw, pixw + 1];
        nbatch = length(idbatch) - 1;
        
        for ib = 1: nbatch
            data = m_in.reg(:, idbatch(ib): idbatch(ib + 1) - 1, :);
            
            switch nflag
                case 1
                    dtmp1 = data .* ~mask(:, idbatch(ib): idbatch(ib + 1) - 1);
                    dtmp1 = dtmp1 .^ 4;
                    dtmp2 = data .* mask(:, idbatch(ib): idbatch(ib + 1) - 1);
                    data = dtmp1 + dtmp2;
                case 2
                    %%% 1. get useful pixels %%%
                    id = mask(:, idbatch(ib): idbatch(ib + 1) - 1);
                    data = reshape(data, [], nf);
                    datat = data(id, :);
                    
                    %%% 2. clean traces %%%
                    tflag = 1;
                    datat = trace_clean(datat, Fs, tflag);
                    
                    %%% put back to raw video %%%
                    data(id, :) = datat;
                    data = reshape(data, pixh, [], nf);
            end
            
            %%% update data %%%
            m_in.reg(:, idbatch(ib): idbatch(ib + 1) - 1, :) = data;
            imax = max(imax, max(max(max(data))));
            imin = min(imin, min(min(min(data))));
        end
        
        ip = 2;
        m_out = normalize_batch(m_in.Properties.Source, 'reg', imax, imin, idbatch, ip);
        
        %%% smooth %%%
        if nflag == 2
            ebatch = ceil(nf / nbatch);
            idbatch = [1: ebatch: nf, nf + 1];
            nbatch = length(idbatch) - 1;
            knl = fspecial('gaussian', round(sqrt(pixh * pixw) / 20), 1);
            for ib = 1: nbatch
                data = m_in.reg(:, :, idbatch(ib): idbatch(ib + 1) - 1);
                data = convn(data, knl, 'same');
                m_in.reg(:, :, idbatch(ib): idbatch(ib + 1) - 1) = data;
            end
        end
    else
        m_out = m_in;
    end
        
    timeh = toc(h);
    disp(['Done noise suppression, time ', num2str(timeh)])
end
