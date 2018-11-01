function Ydcln = dirt_clean(Y, sz, ispara)
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
    Ydcln = zeros(size(Y), class(Y));
    if ispara
        parfor i = 1: nframes
            I = Y(:, :, i);
            tmp = imgaussfilt(I, sz) - I;
            tmp(tmp < 0) = 0;
            Ydcln(:, :, i) = tmp;
            if mod(i, 100) == 0
                disp(['Done frame #', num2str(i), '/', num2str(nframes)])
            end
        end
    else
        for i = 1: nframes
            I = Y(:, :, i);
            tmp = imgaussfilt(I, sz) - I;
            tmp(tmp < 0) = 0;
            Ydcln(:, :, i) = tmp;
            if mod(i, 100) == 0
                disp(['Done #', num2str(i), '/', num2str(nframes)])
            end
        end
    end
end