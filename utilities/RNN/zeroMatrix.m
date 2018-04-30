function M = zeroMatrix(sz, gpu)
%%% Jinghao Lu 02/24/2017 %%%

    if gpu == 1
        M = zeros(sz, 'gpuArray');
    else
        M = zeros(sz);
    end
end