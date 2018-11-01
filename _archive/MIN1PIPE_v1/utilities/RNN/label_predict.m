function predicts = label_predict(P)
%%% Jinghao Lu 02/26/2017 %%%

    [~, predicts] = max(P);
    predicts = predicts' - 1;
end