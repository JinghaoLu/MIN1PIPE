function bg = bg_extract(data)
    %%% get threshold for adaptive background intensity %%%
    [h, w, nf] = size(data);
    n = h * w;
    data = reshape(data, n, nf);
    nuse = 100;
    rg = randsample(n, nuse);
    a = abs(diff(data(rg, :), 1, 2))';
%     thres1 = prctile(a(:), 50);
    [N, hh] = histcounts(a(:));
    [~, idm] = max(N);
    thres1 = hh(idm + 1);

    %%% compute general background intensity %%%
    adif = abs(diff(data, 1, 2));
    tmpt = adif < thres1;
    tmp1 = [tmpt, zeros(n, 1)];
    tmp2 = [zeros(n, 1), tmpt];
    clear tmpt
    idall = tmp1 | tmp2;
    clear tmp1 tmp2
    thres2 = min(data, [], 2) + range(data, 2) / 2;
    idall = idall & (data > thres2);
    bg = sum(data .* idall, 2) ./ sum(idall, 2);
    
    %%% compute NaN pixels %%%
    idnan = find(isnan(bg));
    bg(idnan) = mean(data(idnan, :), 2);
    bg = reshape(bg, h, w);
end

