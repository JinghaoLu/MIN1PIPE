function xbatch = batch_init(xin, bsize, lbin)
%%% Jinghao Lu 02/28/2017 %%%

    nsample = size(xin, 1);
    nbatch = ceil(nsample / bsize);
    lbp = lbin == 1;
    lbn = lbin == 0;
    xpos = xin(lbp, :);
    xneg = xin(lbn, :);
    [~, idxp] = sort(rand(1, size(xpos, 1)));
    [~, idxn] = sort(rand(1, size(xneg, 1)));
    xbatch.xin = cell(1, nbatch);
    xbatch.lbin = cell(1, nbatch);
    for i = 1: nbatch
        rg = [(i - 1) * bsize / 2 + 1, min(nsample / 2, i * bsize / 2)];
        xbatch.xin{i} = [xpos(idxp(rg(1): rg(2)), :); xneg(idxn(rg(1): rg(2)), :)];
        xbatch.lbin{i} = [lbp(idxp(rg(1): rg(2))); lbn(idxn(rg(1): rg(2)))];
    end
end