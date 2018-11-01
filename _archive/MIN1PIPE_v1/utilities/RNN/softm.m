function P = softm(Ut, lstmc)
%%% Jinghao Lu 02/24/2017 %%%

    for l = 1: length(Ut)
        P = exp(Ut{l} * lstmc{end, end}.ht);
        P = bsxfun(@rdivide, P, sum(P, 1));
    end
end