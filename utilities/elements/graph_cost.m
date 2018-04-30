function wts = graph_cost(aint, bint)
% image graph cost function
%   Jinghao Lu, 06/11/2016

    alpha = 400; %%% should be large enough for the range (0, 1) %%%
    wts = exp(alpha * (aint - bint));
%     wts = 0.3 * exp(-alpha * aint) + 0.3 * exp(-alpha * bint) + 0.4 * exp(alpha * abs(aint - bint));
end