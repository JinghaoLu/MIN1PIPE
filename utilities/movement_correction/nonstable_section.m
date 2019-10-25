function m = nonstable_section(m, sttn, stpn, se, pixs, scl, sigma_x, sigma_f, sigma_d, maskc)
% register frames within nonstable sections
%   Jinghao Lu, 02/02/2018

    [pixh, pixw, nf] = size(m, 'reg');
    if nargin < 4 || isempty(se)
        defpar = default_parameters;
        se = defpar.neuron_size;
    end
    
    if nargin < 5 || isempty(pixs)
        pixs = min(pixh, pixw);
    end
    
    if nargin < 6 || isempty(scl)
        defpar = default_parameters;
        scl = defpar.mc_scl;
    end
    
    if nargin < 7 || isempty(sigma_x)
        defpar = default_parameters;
        sigma_x = defpar.mc_sigma_x;
    end
    
    if nargin < 8 || isempty(sigma_f)
        defpar = default_parameters;
        sigma_f = defpar.mc_sigma_f;
    end
    
    if nargin < 9 || isempty(sigma_d)
        defpar = default_parameters;
        sigma_d = defpar.mc_sigma_d;
    end
    
    if nargin < 10 || isempty(maskc)
        maskc = true(pixh, pixw);
    end
    
    %% generate nonstable clusters %%
    %%% compute the range of the nonstable sections %%%
    [bstt, bstp] = nonstable_range(sttn, stpn, nf);
    
    if ~isempty(bstt)
        %% further adjust sections for memory fitness: single nonstable section size %%
        alen = bstp - bstt + 1;
        nff = max(alen);
        ttype = class(m.reg(1, 1, 1));
        stype = parse_type(ttype);
        nsize = pixh * pixw * nff * stype; %%% size of single %%%
        nbatch = batch_compute(nsize);
        ebatch = round(nff / nbatch);
        idx = find(alen > ebatch);
        for i = 1: length(idx)
            tmp = bstt(idx(i)): ebatch: stpn(idx(i));
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
        nsize = pixh * pixw * nff * 8 * feature('numCores') * 4; %%% 8: size of double; 4: rough running size of logdemons %%%
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
        flag = false;
        
        for i = 1: nbatch
            %%% data preparation %%%
            %%%% get the frames of the nonstable sections %%%%
            regpara = cell(1, idbatch(i + 1) - idbatch(i));
            for ii = idbatch(i) + 1: idbatch(i + 1)
                regpara{ii - idbatch(i)} = m.reg(1: pixh, 1: pixw, max(1, bstt(ii) - 1): min(nf, bstp(ii)));
            end
            
            %%%% augment edge blocks %%%%
            if bstt(1) == 1
                if i == 1
                    regpara{1} = cat(3, m.reg(1: pixh, 1: pixw, sttn(1)), regpara{1});
                end
            end
            
            %%% LK-LogDemons %%%
            stof = idbatch(i);
            parfor ii = idbatch(i) + 1: idbatch(i + 1)
                regcur = double(regpara{ii - stof});
                
                %%%% round 1 logdemons %%%%
                regcur = lk_logdemons_unit(regcur, se, pixs, scl, sigma_x, sigma_f, sigma_d, flag, maskc);
                
                %%%% compute similarity scores %%%%
                tmp = reshape(regcur, pixh * pixw, []);
                cscores = diag(corr(tmp), 1);
                idc = find(cscores < 0.5);
                
                for j = 1: length(idc)
                    tmp1 = tmp(:, idc(j));
                    tmp2 = tmp(:, idc(j) + 1);
                    sct1 = (abs(norm_inner(tmp1', tmp1 - tmp2)) + abs(norm_inner(tmp1', tmp2 - tmp1))) / 2;
                    sct2 = (abs(norm_inner(tmp2', tmp1 - tmp2)) + abs(norm_inner(tmp2', tmp2 - tmp1))) / 2;
                    if sct1 > sct2
                        regcur(:, :, idc(j) + 1) = regcur(:, :, idc(j));
                    else
                        regcur(:, :, idc(j)) = regcur(:, :, idc(j) + 1);
                    end
                end
                
                %%%% round 2 logdemons if frame replaced %%%%
                if ~isempty(idc)
                    regcur = lk_logdemons_unit(regcur, se, pixs, scl, sigma_x, sigma_f, sigma_d, flag, maskc);
                end
                
                %%%% update the variable %%%%
                regpara{ii - stof} = regcur;
                
                %%%% release worker memory %%%%
                regcur = [];
                tmp = [];
                
                %%%% stats display %%%%
                if length(bstt) < 10
                    disp(['Done nonstable-LogDemons section # ', num2str(ii + stof), '/', num2str(length(bstt))])
                else
                    if mod(ii, round(length(bstt) / 10)) == 0
                        disp(['Done nonstable-LogDemons section # ', num2str(ii + stof), '/', num2str(length(bstt))])
                    end
                end
                %             disp(['Done nonstable-LogDemons section # ', num2str(ii), '/', num2str(length(bstt))])
            end
            
            %%% combine nonstable sections into the full tensor %%%
            for ii = idbatch(i) + 1: idbatch(i + 1)
                m.reg(1: pixh, 1: pixw, bstt(ii): bstp(ii)) = regpara{ii - idbatch(i)}(:, :, 2: end);
            end
        end
    end
end
