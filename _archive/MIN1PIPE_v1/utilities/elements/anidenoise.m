function YDeN = anidenoise(Y, sz, ispara, iter, dt, kappa, opt)
% Yblur = anidenoise(Y) denoise using anisotropic diffusion
%   denoise image/sequence
%   Y is input image/sequence
%   Yblur is output image/sequence
%   Jinghao Lu 05/17/2016

%     hani = tic;
    %% initialization %%    
    %%% initialize parameters %%%
    if nargin < 2 || isempty(sz)
        defpar = default_parameters;
        sz = defpar.neuron_size;
    end
    
    if nargin < 3 || isempty(ispara)
        ispara = 0;
    end
    
    if nargin < 4 || isempty(iter)
        defpar = default_parameters;
        iter = defpar.anidenoise_iter; 
    end
    
    if nargin < 5 || isempty(dt)
        defpar = default_parameters;
        dt = defpar.anidenoise_dt; %%% maximum for numerical stability, and reduce iterations %%%
    end
    
    if nargin < 6 || isempty(kappa)
        defpar = default_parameters;
        kappa = defpar.anidenoise_kappa; %%% any value above 0.1 for normalized image %%%
    end
    
    if nargin < 7 || isempty(opt)
        defpar = default_parameters;
        opt = defpar.anidenoise_opt;
    end
    
    %%% prepare data %%%
    Y = padarray(Y, [sz, sz], 'replicate');
    nframes = size(Y, 3);
    YDeN = zeros(size(Y), class(Y));
        
%     %%% examine and remove bad pixels %%%
%     temp = max(Y, [], 3) - min(Y, [], 3);
%     [x, y] = find(temp > max(0.2, 20 * std(temp(:))));
%     for i = 1: length(x)
%         Y(x(i), y(i), :) = 0;
%     end
    
    %% begin anisotropic diffusion %%
    if nframes == 1
        YDeN = anisodiff_unit(Y, iter, dt, kappa, opt);   
    else
        if ispara
            if isempty(gcp('nocreate'))
                parpool(feature('numCores'));
            end
            parfor i = 1: nframes
                I = Y(:, :, i);
                YDeN(:, :, i) = anisodiff_unit(I, iter, dt, kappa, opt);
                if mod(i, 100) == 0
                    disp(['done frame #', num2str(i), '/', num2str(nframes)])
                end
            end
        else
            for i = 1: nframes
                I = Y(:, :, i);
                YDeN(:, :, i) = anisodiff_unit(I, iter, dt, kappa, opt);
                if mod(i, 100) == 0
                    disp(['done frame #', num2str(i), '/', num2str(nframes)])
                end
            end
        end
    end
    
%     time = toc(hani);
%     disp(['Done anidenoise ', num2str(time), ' seconds'])    
end