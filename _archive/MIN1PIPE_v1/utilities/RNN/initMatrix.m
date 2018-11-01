function out = initMatrix(h, w, gpu)
%%% Jinghao Lu 02/25/2017 %%%

    amp = 0.1;
    if gpu == 1
        out = amp * (rand([h, w], 'double', 'gpuArray') - 0.5);
    else
        out = amp * (rand([h, w]) - 0.5);
    end
end