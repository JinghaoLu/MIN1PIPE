function lstmc = fprop(xin, Wt, ht1, ct1)
%%% Jinghao Lu 02/26/2017 %%%

    normp = [1, 0];
    [bsize, tlen] = size(xin);
    layern = length(Wt);
    hdim4 = size(Wt{1}, 1);
    hdim = hdim4 / 4;
    Wi = cell(1, layern);
    Wf = cell(1, layern);
    Wo = cell(1, layern);
    Wc = cell(1, layern);
    lstmc = cell(layern, tlen);
    for l = 1: layern
        Wi{l} = Wt{l}(1: hdim, :);
        Wf{l} = Wt{l}(hdim + 1: 2 * hdim, :);
        Wo{l} = Wt{l}(2 * hdim + 1: 3 * hdim, :);
        Wc{l} = Wt{l}(3 * hdim + 1: end, :);
        ht1{l} = ht1{l}(:, 1: bsize);
        ct1{l} = ct1{l}(:, 1: bsize);
    end
    
    for i = 1: tlen
        datain = xin(:, i)';
        for l = 1: layern
            lstmc{l, i}.incat = [datain; ht1{l}; ones(1, bsize)]; %%% x + h + b %%%
            lstmc{l, i}.igate = sigmf(Wi{l} * lstmc{l, i}.incat, normp);
            lstmc{l, i}.fgate = sigmf(Wf{l} * lstmc{l, i}.incat, normp);
            lstmc{l, i}.ogate = sigmf(Wo{l} * lstmc{l, i}.incat, normp);
            lstmc{l, i}.c_hat = tanh(Wc{l} * lstmc{l, i}.incat);
            lstmc{l, i}.ct = lstmc{l, i}.fgate .* ct1{l} + lstmc{l, i}.igate .* lstmc{l, i}.c_hat;
            lstmc{l, i}.tct = tanh(lstmc{l, i}.ct);
            lstmc{l, i}.ht = lstmc{l, i}.ogate .* lstmc{l, i}.tct;
            ct1{l} = lstmc{l, i}.ct;
            ht1{l} = lstmc{l, i}.ht;
            datain = ht1{l};
        end
    end
end