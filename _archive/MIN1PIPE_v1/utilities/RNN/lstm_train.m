function [Wall, Uall, tp, acc, cren, acct] = lstm_train(xin, lbin, pars, xval, lbval, parsval, optim, Wall, Uall)
%%% Jinghao Lu 03/02/2017 %%%

    %% initialize %%
    if nargin < 7
        optim = 'adam';
    end
    optimizer = optimizer_init(pars, optim);
    nsample = length(lbin);
    [Wt, Ut, ht1, ct1] = lstm_init(pars);
    if nargin >= 8
        Wt = Wall{end};
        Ut = Uall{end};
    end
    Wall = {};
    Uall = {};
    [~, ~, ht1val, ct1val] = lstm_init(parsval);
    tp = [];
    cren = [];
    acc = [];
    acct = [];
    flag = 0;
    T = 0;
    lr = 0;
    
    %% minibatch-based %%
    nepoch = pars.nepoch;
    nbatch = ceil(nsample / pars.bsize);
    
    for ie = 1: nepoch
        xbatch = batch_init(xin, pars.bsize, lbin);
        for ib = 1: nbatch
            %%% forward prop %%%
            lstmc = fprop(xbatch.xin{ib}, Wt, ht1, ct1);
            
            %%% softmax %%%
            P = softm(Ut, lstmc);
            acct = [acct, gather(lstm_val(xbatch.lbin{ib}, P))];
            
            %%% cross-entropy loss %%%
            crencost = cross_entropy_loss(P, xbatch.lbin{ib});
            tp = [tp, crencost];
            disp(['Cross Entropy Loss ', num2str(crencost), ', acc ', num2str(acct(end)), ', ---- epoch ', num2str(ie), ', iter ', num2str(ib)])

            %%% BPTT %%%
            gradt = bprop(lstmc, Wt, Ut, P, xbatch.lbin{ib}, pars);
            
            %%% check grad %%%
            check_grad(gradt.dLdW, Wt, Ut, xbatch.xin{ib}, xbatch.lbin{ib}, ht1, ct1, pars.check)

            %%% SGD %%%
            t = (ie - 1) * nbatch + ib;
            if strcmp(optim, 'adam')
%                 lr = 0.001 * (1 / (1 + 1 * ie));
%                 lr = 0.001 / (1 + 0.2 * t);
%                 lr = 0.001 * exp(-0.06 * t);
                if flag == 1 && length(cren) > 5 && mean(abs(diff(cren(max(1, end - 4): end)))) < 0.0001
                    T = t;
                    pars.lrate1 = pars.lrate2;
                    flag = 0;
                else
                    lr = pars.lrate1(t - T, lr);
                end
                [Wt, optimizer.Wm, optimizer.Wv] = sgd_adam(Wt, gradt.dLdW, optimizer.Wm, optimizer.Wv, t, lr); %%% t is sgd step %%%
                [Ut, optimizer.Um, optimizer.Uv] = sgd_adam(Ut, gradt.dLdU, optimizer.Um, optimizer.Uv, t, lr);
            elseif strcmp(optim, 'rmsprop')
                [Wt, optimizer.W] = sgd_rmsprop(Wt, gradt.dLdW, optimizer.W); 
                [Ut, optimizer.U] = sgd_rmsprop(Ut, gradt.dLdU, optimizer.U); 
            end
%             disp(num2str(optimizer{1}(100,100)))
            
            for l = 1: pars.layern
                Wall{t}{l} = gather(Wt{l});
            end
            Uall{t}{1} = gather(Ut{1});
            
            %%% validation %%%
            lstmc = fprop(xval, Wt, ht1val, ct1val);
            P = softm(Ut, lstmc);
            acc = [acc, gather(lstm_val(lbval, P))];
            crencost = cross_entropy_loss(P, lbval);
            cren = [cren, crencost];
            disp(['acc ', num2str(acc(end)), ', cross entropy ', num2str(crencost)])
            if acc(end) > 1
                return
            end
        end
    end
end









