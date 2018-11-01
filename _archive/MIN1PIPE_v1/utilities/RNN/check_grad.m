function check_grad(grad, W, U, xin, lbin, ht1, ct1, ifcheck)
%%% Jinghao Lu 02/25/2017 %%%

    if ifcheck
        e = 1e-4;
        i = 100;
        j = 100;
        l = 1;
        grad1 = grad{l}(i, j);
        W{l}(i, j) = W{l}(i, j) + e;
        lstm = fprop(xin, W, ht1, ct1);
        PP = softm(U, lstm);
        cost1 = cross_entropy_loss(PP, lbin);
        
        W{l}(i, j) = W{l}(i, j) - 2 * e;
        lstm = fprop(xin, W, ht1, ct1);
        PP = softm(U, lstm);
        cost2 = cross_entropy_loss(PP, lbin);
        
        grad2 = (cost1 - cost2) / (2 * e);
        disp([num2str(grad1), ', ', num2str(grad2), ', ', num2str(grad1 - grad2)])
    end
end