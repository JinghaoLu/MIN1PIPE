function [labelout, prob] = seeds_cleansing_rnn(xtest, model, stride)
    
    %%% for entire time seris %%%
    if nargin < 3
        stride = 3;
    end
    [nsamples, seqlen] = size(xtest);
    tstep = ceil((seqlen - model.tlen + 1) / stride);
    tlen = model.tlen;
    gpu = model.pars.gpu;
    pars = model.pars;
    nhid = size(model.W{1}, 1) / 4;
    if gpu
        W = cell(1, pars.layern);
        for l = 1: pars.layern
            W{l} = gpuArray(model.W{l});
        end
        U = {gpuArray(model.U{1})};
    else
        W = cell(1, pars.layern);
        for l = 1: pars.layern
            W{l} = gather(model.W{l});
        end
        U = {gather(model.U{1})};
    end
   
    gpucurr = gpuDevice;
    tltmem = gpucurr.AvailableMemory;
    cbatch = max(1, floor(tltmem / (2 * 70 * nsamples * nhid * tlen)));
    tclust = ceil(tstep / cbatch);
    labelout = zeros(nsamples, tstep);
    prob = zeros(2, nsamples, tstep);
    for ic = 1: tclust
        disp(['Begin tcluster ', num2str(ic), '/', num2str(tclust)])
        trange = [(ic - 1) * cbatch + 1, min(ic * cbatch, tstep)];
        tstepcur = diff(trange) + 1;
        data = zeroMatrix([nsamples, tlen, tstepcur], gpu);
        for j = trange(1): trange(2)
            strbg = (j - 1) * stride;
            data(:, :, j - trange(1) + 1) = xtest(:, strbg + 1: strbg + tlen);
            data(:, :, j - trange(1) + 1) = bsxfun(@minus, data(:, :, j - trange(1) + 1), min(data(:, :, j - trange(1) + 1), [], 2));
        end
        data = permute(data, [3, 1, 2]);
        data = reshape(data, tstepcur * nsamples, tlen);

        gpucurr = gpuDevice;
        tltmem = gpucurr.AvailableMemory;
        batchstp = floor(tltmem / (2 * 70 * nhid * tlen));
        pars.bsize = batchstp;
        nbatch = ceil(size(data, 1) / batchstp);
        labeloutt = cell(nbatch, 1);
        probt = cell(1, nbatch);
        xin = cell(1, nbatch);
        for i = 1: nbatch
            xin{i} = data((i - 1) * batchstp + 1: min(tstepcur * nsamples, i * batchstp), :);
            xin{i} = bsxfun(@minus, xin{i}, min(xin{i}, [], 2));
        end
        data = xin;
        [~, ~, ht1test, ct1test] = lstm_init(pars);
        for i = 1: nbatch
            [labeloutt{i}, probt{i}] = lstm_predict_unit(data{i}, W, U, ht1test, ct1test);
            disp(['----- tcluster ', num2str(ic), ', batch #', num2str(i), '/', num2str(nbatch)])
        end
        labeloutt = cell2mat(labeloutt);
        labeloutt = reshape(labeloutt, tstepcur, nsamples);
        labelout(:, trange(1): trange(2)) = labeloutt';

        probt = cell2mat(probt);
        probt = reshape(probt, 2, [], nsamples);
        prob(:, :, trange(1): trange(2)) = permute(probt, [1, 3, 2]);
    end
    
    labelout = max(labelout, [], 2);
    prob = squeeze(max(prob, [], 3));
end