function Ydebg = bg_remove(Y, sz, ispara)
% Ydebg = bg_remove(Y, ispara, sz) remove background using morphological 
% opening
%   Jinghao Lu, 01/16/2016

%     hbg = tic;
    %% initialization %%
    %%% initialize parameters %%%
    if nargin < 2 ||isempty(sz)
        defpar = default_parameters;
        sz = defpar.neuron_size;
    end
    
    if nargin < 3 || isempty(ispara)
        ispara = 0;
    end
    
    %%% prepare parallel computing %%%
    if ispara
        if isempty(gcp('nocreate'))
            parpool(feature('numCores'));
        end
    end
    
    %% begin removing background %%
    nframes = size(Y, 3);
    Ydebg = zeros(size(Y(sz + 1: end - sz, sz + 1: end - sz, :)), class(Y));
    k = strel('disk', sz);
    if ispara
        parfor i = 1: nframes
            I = Y(:, :, i);
            bg = imopen(I, k);
            tmp = I - bg;
%             tmp(1: sz, :) = 0;
%             tmp(:, end - sz: end) = 0;
%             tmp(end - sz: end, :) = 0;
%             tmp(:, 1: sz) = 0;
            Ydebg(:, :, i) = tmp(sz + 1: end - sz, sz + 1: end - sz);
            if mod(i, 100) == 0
                disp(['Done frame #', num2str(i), '/', num2str(nframes)])
            end
        end
    else
        for i = 1: nframes
            I = Y(:, :, i);
            bg = imopen(I, k);
            tmp = I - bg;
%             tmp(1: sz, :) = 0;
%             tmp(:, end - sz: end) = 0;
%             tmp(end - sz: end, :) = 0;
%             tmp(:, 1: sz) = 0;
            Ydebg(:, :, i) = tmp(sz + 1: end - sz, sz + 1: end - sz);
            if mod(i, 100) == 0
                disp(['Done #', num2str(i), '/', num2str(nframes)])
            end
        end
    end
    
%     %%% norm after remove %%%
%     Ydebg = normalize(Ydebg);
%     time = toc(hbg);
%     disp(['Done bg remove, total time: ', num2str(time), ' seconds'])
end