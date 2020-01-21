function [gx, gy] = gradient_fast(I)
% inlining the gradient computing
%   Jinghao Lu, 05/02/2019

    gx = I;
    gy = gx;
    gy(1, :) = I(2, :) - I(1, :);
    gy(end, :) = I(end, :) - I(end - 1, :);
    gy(2: end - 1, :) = (I(3: end, :) - I(1: end - 2, :)) / 2;
    gx(:, 1) = I(:, 2) - I(:, 1);
    gx(:, end) = I(:, end) - I(:, end - 1);
    gx(:, 2: end - 1) = (I(:, 3: end) - I(:, 1: end - 2)) / 2;
end
