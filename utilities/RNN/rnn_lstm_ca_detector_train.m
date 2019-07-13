pars.layern = 2;
pars.indim = 1;
pars.hdim = 200;
pars.outdim = 2;
pars.bsize = 200;
pars.nepoch = 100;
pars.gpu = 1;
pars.check = 0;
pars.lrate1 = @(t, lr) 0.001 * exp(-0.06 * t);
% pars.lrate1 = @(t, lr) 0.001 ./ (1 + 0.2 * t); 
pars.lrate2 = @(t, lr) lr ./ (1 + 0.02 * t); 

% idx = find(train_lbl == 1);
% npos = length(idx);
% npostrain = min(512, npos * 0.9);
% xin = train_seq(idx(1: npostrain), :);
% lbin = train_lbl(idx(1: npostrain), :);
% xval = train_seq(idx(npostrain + 1: end), :);
% lbval = train_lbl(idx(npostrain + 1: end));
% idx = find(train_lbl == 0);
% xin = [xin; train_seq(idx(1: npostrain), :)];
% lbin = [lbin; train_lbl(idx(1: npostrain), :)];
% xval = [xval; train_seq(idx(npostrain + 1: npos), :)];
% lbval = [lbval; train_lbl(idx(npostrain + 1: npos))];
xin = train_seq;
lbin = train_lbl;
xval = train_seq_val;
lbval = train_lbl_val;
if pars.gpu == 1
    gpuDevice(1);
    xin = gpuArray(xin);
    lbin = gpuArray(lbin);
    xval = gpuArray(xval);
    lbval = gpuArray(lbval);
end

parsval = pars;
parsval.bsize = length(lbval);

tic
% [model.W, model.U, tp, acc, cren, acct] = lstm_train(xin, lbin, pars, xval, lbval, parsval, 'adam', model.W, model.U);
% [model.W, model.U, tp, acc, cren, acct] = lstm_train(xin, lbin, pars, xval, lbval, parsval, 'rmsprop');
[model.W, model.U, tp, acc, cren, acct] = lstm_train(xin, lbin, pars, xval, lbval, parsval, 'adam');
toc


%%% complete model %%%
modelbu = model;
model.W = modelbu.W{end};
model.U = modelbu.U{end};
model.tlen = size(xin, 2);
model.pars = pars;
% seq = train_seq;
% lbl = train_lbl;
seq = sig;
% lbl = test_lbl;
[labelout, P] = seeds_cleansing_rnn(seq, model);

PP = reshape(P, 2, [], size(seq, 1));
aa = squeeze(max(PP, [], 2));
lblb = sum(labelout, 1) > 0;
tpr = sum(lblb & lbl') / sum(lbl);
fpr = sum(lblb &(~lbl')) / (length(lblb) - sum(lbl));


