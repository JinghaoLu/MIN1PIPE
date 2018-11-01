function acc = lstm_val(lbval, P)
%%% Jinghao Lu 02/26/2017 %%%

    predicts = label_predict(P);
    acc = sum((predicts - lbval) == 0) / length(lbval);
end