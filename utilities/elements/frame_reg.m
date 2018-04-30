function [regfn, acorrf, acorr, scl] = frame_reg(Y, Fs, pixs, scl, sigma_x, sigma_f, sigma_d)
% register movies with the hierarchical movement correction
%   Jinghao Lu, 09/01/2017

    hreg = tic;
    %% initialization %%
    %%% initialize parameters %%%
    if nargin < 2 || isempty(Fs)
        defpar = default_parameters;
        Fs = defpar.Fsi_new;
    end
    
    if nargin < 3 || isempty(pixs)
        [pixh, pixw, ~] = size(Y);
        pixs = min(pixh, pixw);
    end
    
    if nargin < 4 || isempty(scl)
        defpar = default_parameters;
        scl = defpar.mc_scl;
    end
    
    if nargin < 5 || isempty(sigma_x)
        defpar = default_parameters;
        sigma_x = defpar.mc_sigma_x;
    end
    
    if nargin < 6 || isempty(sigma_f)
        defpar = default_parameters;
        sigma_f = defpar.mc_sigma_f;
    end
    
    if nargin < 7 || isempty(sigma_d)
        defpar = default_parameters;
        sigma_d = defpar.mc_sigma_d;
    end

    %%% prepare parallel computing %%%
    if isempty(gcp('nocreate'))
        parpool(feature('numCores'));
    end
        
    %%% preprocess Y first %%%
    dthres = 0.1;
    mskpre = dominant_patch(Y, dthres);
    Y = Y .* mskpre;

    %%% get translation score %%%
    fprintf('Begin initial computation of translation score \n')
    acorr = get_trans_score(Y, [], [], 1);

    %%% cluster movie into hierarchical stable-nonstable sections %%%
    [stt, stp, flag, scl] = hier_clust(acorr, Fs, pixs, scl); %%% flag: real or fake clusters %%%
    time = toc(hreg);
    fprintf(['Done initialization, ', num2str(time), ' seconds \n'])

    %% intra-section registration %%
    fprintf('Begin intra-section \n')
    reg_intra = intra_section(Y, stt, stp, pixs, scl, sigma_x, sigma_f, sigma_d, flag);
    time = toc(hreg);
    fprintf(['Done intra-section, ', num2str(time), ' seconds \n'])
    clear Y

    %% inter-section registration %%
    fprintf('Begin inter-section ... ')
    [reg_inter, ~, ~] = inter_section(reg_intra, stt, stp, pixs, scl, sigma_x, sigma_f, sigma_d);
    time = toc(hreg);
    fprintf(['Done inter-section, ', num2str(time), ' seconds \n'])
    clear reg_intra

    %% nonstable-section registration %%
    fprintf('Begin nonstable-section \n')
    reg_nonstable = nonstable_section(reg_inter, stt, stp, pixs, scl, sigma_x, sigma_f, sigma_d);
    time = toc(hreg);
    fprintf(['Done nonstable-section, ', num2str(time), ' seconds \n'])
    clear reg_inter
        
    %% final preparation for output %%
    %%% final score %%%
    disp('Begin final computation of translation score')
    acorrf = get_trans_score(reg_nonstable, [], [], 1);
    regfn = reg_nonstable;
    time = toc(hreg);
    disp(['Done frame reg, total time: ', num2str(time), ' seconds'])
end

