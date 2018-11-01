function [Wt, m, v] = sgd_adam(Wt1, dW, mt1, vt1, t, lr, beta)
%%% Jinghao Lu 02/28/2017 %%%

    if nargin < 6
        lr = 0.001;
    end
    
    if nargin < 7
        beta = [0.9, 0.999];
    end
    
    layern = length(dW);
    m = cell(1, layern);
    v = cell(1, layern);
    Wt = cell(1, layern);
    epsilon = eps;
    beta1 = beta(1);
    beta2 = beta(2);
    
    for l = 1: layern
        %%% Update biased moment estimate %%%
        m{l} = beta1 .* mt1{l} + (1 - beta1) .* dW{l};
        v{l} = beta2 .* vt1{l} + (1 - beta2) .* (dW{l} .^ 2);
        
        %%% Compute bias-corrected moment estimate %%%
        m_hat = m{l} ./ (1 - beta1 ^ t);
        v_hat = v{l} ./ (1 - beta2 ^ t);
        
        %%% Update decision variables %%%
        Wt{l} = Wt1{l} - lr .* m_hat ./ (sqrt(v_hat) + epsilon);
    end
end