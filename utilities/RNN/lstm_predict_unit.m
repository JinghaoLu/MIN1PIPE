function [labelout, P] = lstm_predict_unit(data, W, U, ht1test, ct1test)
%%% Jinghao Lu 03/03/2017 %%%

    lstmc = fprop(data, W, ht1test, ct1test);
    P = gather(softm(U, lstmc));
    labelout = gather(label_predict(P));
%     lbpred = reshape(lbpred, tstep, size(temp, 1) / tstep);
%     labelout = sum(lbpred, 1) > 0;
%     labelout = lbpred;
end