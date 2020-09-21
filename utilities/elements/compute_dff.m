function dff = compute_dff(sigfn, bgfn, bgffn, seedsfn)
    mn = zeros(size(sigfn, 1), 1);
    for i = 1: size(sigfn, 1)
        tmp = sigfn(i, :);
        edges = linspace(min(tmp), max(tmp), 101);
        n = ksdensity(tmp, edges);
        [~, id] = max(n);
        mn(i) = edges(id);
    end
    dff = sigfn ./ (mean(bgffn) * bgfn(seedsfn) + mn);
end