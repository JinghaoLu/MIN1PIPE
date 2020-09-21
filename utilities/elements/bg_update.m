function [b, bf] = bg_update(m, roi, sig)
    %% background update %%
    [d1, d2, d3] = size(m, 'reg');
    d = d1 * d2;
    nsize = d * d3 * 8; %%% size of double %%%
    nbatch = batch_compute(nsize);
    ebatch = ceil(d / nbatch);
    eb = floor(ebatch / d1);
    idbatch = [1: eb: d2, d2 + 1];
    recbatch = (idbatch - 1) * d1 + 1;
    nbatch = length(idbatch) - 1;
    
    bf = zeros(1, d3);
    denom = 0;
    for i = 1: nbatch
        tmp = m.reg(1: d1, idbatch(i): idbatch(i + 1) - 1, 1: d3);
        tmp = reshape(tmp, [], d3);
        tmpp = roi(recbatch(i): recbatch(i + 1) - 1, :) * sig;
        mn = mean(tmp, 2) - mean(tmpp, 2);
        mask = mn > 0;
        tmp = tmp(mask, :);
        tmpp = tmpp(mask, :);
        bf = bf + sum(tmp, 1) - sum(tmpp, 1);
        denom = denom + sum(mask);
    end
    bf = max(0, double(bf) / denom);
    
    b = zeros(d, 1);
    for i = 1: nbatch
        tmp = m.reg(1: d1, idbatch(i): idbatch(i + 1) - 1, 1: d3);
        tmp = double(reshape(tmp, d1 * (idbatch(i + 1) - idbatch(i)), d3));
        Yf = tmp * bf';
        tmpp = roi(recbatch(i): recbatch(i + 1) - 1, :) * sig;
        b(d1 * (idbatch(i) - 1) + 1: d1 * (idbatch(i + 1) - 1)) = max((Yf - tmpp * bf') / (bf * bf'), 0);
    end
end