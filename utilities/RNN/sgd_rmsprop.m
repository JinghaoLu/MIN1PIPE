function [Wt, et] = sgd_rmsprop(Wt1, dW, et1, lr, gamma)
%%% Jinghao Lu 03/02/2017 %%%

    if nargin < 4
        lr = 0.0001;
    end
    
    if nargin < 5
        gamma = 0.9;
    end
    
    layern = length(dW);
    et = cell(1, layern);
    Wt = cell(1, layern);
    epsilon = eps;
    
    for l = 1: layern
        %%% Update past average estimate %%%
        et{l} = gamma .* et1{l} + (1 - gamma) .* (dW{l} .^ 2);
                
        %%% Update decision variables %%%
        Wt{l} = Wt1{l} - lr .* dW{l} ./ (sqrt(et{l} + epsilon));
    end
end