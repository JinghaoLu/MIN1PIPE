function [roi, sig, seeds, datasmth, cutoff, pkcutoff] = final_seeds_select(m, roi, sig, seeds, datasmth, cutoff, pkcutoff, sz, maxall)
    [pixh, pixw, nf] = size(m, 'reg');    

    %% kstest %%
    id1 = false(1, length(seeds));
    for i = 1: length(seeds)
        id1(i) = kstest(sig(i, :));
    end
    
    %% skewness %%
    skw = skewness(sig');
    id2 = skw > 0;
    
    %% ksdensity distribution %%
    id3 = false(1, length(seeds));
    id4 = false(1, length(seeds));
    cthres = 60;
    ithres = 60;
    kss = zeros(length(seeds), 100);
    parfor i = 1: length(seeds)
        rg = linspace(min(sig(i, :)), max(sig(i, :)), 100);
        tsig = sig(i, :);
        tsig = tsig(tsig > prctile(tsig, 1) & tsig < prctile(tsig, 99));
        tmp = ksdensity(tsig, rg);
        [~, idmax] = max(tmp);
        id3(i) = idmax < ithres;
        tmp = cumsum(tmp) / sum(tmp);
        kss(i, :) = tmp;
        tmp = find(tmp > 0.5, 1);
        id4(i) = tmp < cthres;
    end
    
    %% intensity filter %%
    ithres = intensity_filter(maxall) * 0.5;
    mtmp = maxall(seeds);
    id5 = mtmp > ithres;
    
    %% spatial filter %%
    krnt = fspecial('gaussian', [pixh, pixw], floor(sz / 2));
    krn = krnt(floor(pixh / 2 - sz): ceil(pixh / 2 + sz), floor(pixw / 2 - sz): ceil(pixw / 2 + sz));
    maxcorr = normxcorr2(krn, maxall);
    [hk, wk] = size(krn);
    maxcorr = maxcorr(ceil(hk / 2): end - ceil(hk / 2) + 1, ceil(wk / 2): end - ceil(wk / 2) + 1);
    [~, id6, ~] = intersect(seeds, find(maxcorr(:) > 0));

    %% get refined variables %%
    id = intersect(find(id1 & id2 & id3 & id4 & id5), id6);
    roi = roi(:, id);
    sig = sig(id, :);
    seeds = seeds(id);
    datasmth = datasmth(id, :);
    cutoff = cutoff(id);
    pkcutoff = pkcutoff(id);
    
end