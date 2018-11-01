function pxluse = vld_prd_slct(data, c, pkc, pxlthres)
% pxluse = vld_prd_slct(data, c, pkc, pxlthres) select valid period for time
%   series correlation
%   Jinghao Lu 08/24/2016
    
    %%% initialization %%%
    if nargin < 4
        pxlthres = 100;
    end
    [szh, nf] = size(data);
    if szh == 1
        data = data';
    end
    
    %%% find large calcium spikes %%%
    [pmx, id, ~, ~] = findpeaks(data);
    pkmxvld = find(pmx > pkc);
    idtmp = id(pkmxvld);
    pmxtmp = pmx(pkmxvld);
    [~, idst] = sort(pmxtmp, 'descend');
    idpks = idtmp(idst);
    pxluse = false(1, nf);

    %%% get useful periods for correlation %%%
    if ~isempty(idpks)
        for j = 1: length(idpks)
            tmp1 = data(1: idpks(j));
            tmp2 = data(idpks(j): end);
            id1 = find(diff(tmp1) < 0, 1, 'last');
            if isempty(id1)
                id1 = 1;
            end
            id2 = find(diff(tmp2) > 0, 1, 'first') + idpks(j) - 1;
            if isempty(id2)
                id2 = length(data);
            end
            pxluse(id1: id2) = true;
        end
    else
        [~, pxltmp] = max(data);
        if pxltmp - pxlthres / 2 < 1
            pxluse(1: pxlthres + 1) = true;
        elseif pxltmp + pxlthres / 2 > length(data)
            pxluse(length(data) - pxlthres - 1: length(data)) = true;
        else
            pxluse(pxltmp - pxlthres / 2: pxltmp + pxlthres / 2) = true;
        end
    end

    %%% judge if got enough periods %%%
    if sum(pxluse) < pxlthres
        pkr = setdiff(1: length(pmx), pkmxvld);
        if ~isempty(pkr)
            pmxr = pmx(pkr);
            idr = id(pkr);
            [pmxrs, idrst] = sort(pmxr, 'descend');
            idrs = idr(idrst);
            k = 1;
            while sum(pxluse) < pxlthres && k <= length(idrs)
                tmp1 = data(1: idrs(k));
                tmp2 = data(idrs(k): end);
                id1 = find(diff(tmp1) < 0, 1, 'last');
                id2 = find(diff(tmp2) > 0, 1, 'first') + idrs(k) - 1;
                if isempty(id2)
                    id2 = length(data);
                end
                pxluse(id1: id2) = true;
                k = k + 1;
            end
        else
            
        end
    end
end