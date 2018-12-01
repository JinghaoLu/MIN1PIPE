function [ref, nfref] = ref_select(Y, nfcomp)
% select reference based on image stacks salient seeds pixels
%   Jinghao Lu, 05/11/2017

    %%% find all high intensity "seeds" pixels %%%
    [pixh, pixw, nf] = size(Y);  
    ns = 10;
    basecur = max(Y, [], 3);
    rmax = imregionalmax(basecur);
    ithres = prctile(basecur(:), 98);
    tmp1 = rmax & (basecur > ithres);
    inds = tmp1(:);
    intens = basecur(inds);
    [~, ids] = sort(intens, 'descend');
    tmp = reshape(bsxfun(@times, Y, rmax & (basecur > ithres)), pixh * pixw, nf);
    tmp2 = tmp(inds, :);
    tmp = tmp2(ids(1: min(ns, size(tmp2, 1))), :);
    
    %%% frame number threshold for individual "seeds" %%%
    if nargin < 2
        nfthres = 10;
    else
        nfthres = ceil(nfcomp / size(tmp, 1));
    end
    
    %%% get feature map; individual features are computed from similar
    %%% frames %%%
    nfref = 0;
    mxs = zeros(pixh, pixw, size(tmp, 1));
    for i = 1: size(tmp, 1)
        ta = tmp(i, :);
        [tat, taid] = sort(ta, 'descend');
        thre = 0.5 * (max(tat) + min(tat));
        idt = tat >= thre;
        if sum(idt) > nfthres
            idt = false(1, length(idt));
            idt(1: nfthres) = true;
        end
        nfref = nfref + sum(idt);
        mxs(:, :, i) = sum(Y(:, :, taid(idt)), 3);
    end
    ref = sum(mxs, 3) ./ nfref;
end
