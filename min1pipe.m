function [file_name_to_save, file_name_reg] = min1pipe(Fsi, Fsi_new, spatialr, se, ismc, flag, isvis)
% main_processing
%   need to decide whether to use parallel computing
%   Fsi: raw sampling rate
%   Fsi_new: in use sampling rate
%   spatialr: spatial downsampling factor
%   Jinghao Lu 06/10/2016

    %% configure paths %%
    min1pipe_init;
    
    %% initialize parameters %%
    if nargin < 1 || isempty(Fsi)
        defpar = default_parameters;
        Fsi = defpar.Fsi;
    end
    
    if nargin < 2 || isempty(Fsi_new)
        defpar = default_parameters;
        Fsi_new = defpar.Fsi_new;
    end
    
    if nargin < 3 || isempty(spatialr)
        defpar = default_parameters;
        spatialr = defpar.spatialr;
    end
    
    if nargin < 4 || isempty(se)
        defpar = default_parameters;
        se = defpar.neuron_size;
    end
    
    if nargin < 5 || isempty(ismc)
        ismc = true;
    end
    
    if nargin < 6 || isempty(flag)
        flag = 1;
    end
    
    if nargin < 7 || isempty(isvis)
        isvis = false;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%% parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% user defined parameters %%%                                     %%%
    Params.Fsi = Fsi;                                                   %%%
    Params.Fsi_new = Fsi_new;                                           %%%
    Params.spatialr = spatialr;                                         %%%
    Params.neuron_size = se; %%% half neuron size; 9 for Inscopix and 5  %%%
                            %%% for UCLA, with 0.5 spatialr separately  %%%
                                                                        %%%
    %%% fixed parameters (change not recommanded) %%%                   %%%
    Params.data_cat_ispara = 0;          %%% if parallel (data collect) %%%
    Params.anidenoise_iter = 4;                   %%% denoise iteration %%%
    Params.anidenoise_dt = 1/7;                  %%% denoise step size %%%
    Params.anidenoise_kappa = 0.5;       %%% denoise gradient threshold %%%
    Params.anidenoise_opt = 1;                %%% denoise kernel choice %%%
    Params.anidenoise_ispara = 1;             %%% if parallel (denoise) %%%   
    Params.bg_remove_ispara = 1;    %%% if parallel (backgrond removal) %%%
    Params.mc_scl = 0.004;      %%% movement correction threshold scale %%%
    Params.mc_sigma_x = 3;  %%% movement correction spatial uncertainty %%%
    Params.mc_sigma_f = 0.5;   %%% movement correction fluid reg weight %%%
    Params.mc_sigma_d = 10;%%% movement correction diffusion reg weight %%%
    Params.pix_select_sigthres = 1.5;     %%% seeds select signal level %%%
    Params.pix_select_corrthres = 0.6; %%% merge correlation threshold1 %%%
    Params.refine_roi_ispara = 1;          %%% if parallel (refine roi) %%%
    Params.merge_roi_corrthres = 0.9;  %%% merge correlation threshold2 %%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% get dataset info %%
    [path_name, file_base, file_fmt] = data_info;
    
    hpipe = tic;
    for i = 1: length(file_base)
        %% data cat %%
        ispara = Params.data_cat_ispara; %%% set to 0 if no need of parallel %%%
        Fsi = Params.Fsi;
        Fsi_new = Params.Fsi_new;
        spatialr = Params.spatialr;
        [frame_all, pixh, pixw, nf] = data_cat(path_name, file_base{i}, file_fmt{i}, Fsi, Fsi_new, spatialr, ispara);
%         disp(['Donenf #: ', num2str(nf)])

        %% norm %%
        norm_frame = normalize(frame_all);
        imaxn = max(norm_frame, [], 3);
        if isvis
            data.norm_frame = norm_frame;
        end
        clear frame_all

        %% get rough roi domain %%
        mask = roi_domain(norm_frame);

        %% anisotropic diffusion %%
        ispara = Params.anidenoise_ispara;
        sz = Params.neuron_size;
        iter = Params.anidenoise_iter;
        dt = Params.anidenoise_dt;
        kappa = Params.anidenoise_kappa;
        opt = Params.anidenoise_opt;
        Yblur = anidenoise(norm_frame, sz, ispara, iter, dt, kappa, opt);
        if isvis
            data.Yblur = Yblur;
        end
        clear norm_frame

        %% background remove %%
        ispara = Params.bg_remove_ispara;
        sz = Params.neuron_size;
        Ydebg = bg_remove(Yblur, sz, ispara); %%% strel 9 for 540 * 720 (inscopix); 5 for ucla scope %%%
        imaxy = max(Ydebg, [], 3);
        if isvis
            data.Ydebg = Ydebg;
        end
        clear Yblur

        %% frame register %%
        if ismc
            pixs = min(pixh, pixw);
            Params.mc_pixs = pixs;
            Fsi_new = Params.Fsi_new;
            scl = Params.mc_scl;
            sigma_x = Params.mc_sigma_x;
            sigma_f = Params.mc_sigma_f;
            sigma_d = Params.mc_sigma_d;
            [reg, corr_score, raw_score, scl] = frame_reg(Ydebg, Fsi_new, pixs, scl, sigma_x, sigma_f, sigma_d);
            Params.mc_scl = scl; %%% update latest scl %%%
        else
            reg = Ydebg;
        end
        
        if isvis
            data.reg = reg;
        end
        clear Ydebg
        time1 = toc(hpipe);
        
        %% select pixel %%
        if flag == 1
            sz = Params.neuron_size;
            Fsi_new = Params.Fsi_new;
            sigthres = Params.pix_select_sigthres;
            corrthres = Params.pix_select_corrthres;
            [roi, sig, bg, bgf, seeds, datasmth0, cutoff0, pkcutoff0] = pix_select(reg, mask, sz, Fsi_new, sigthres, corrthres);
            time1 = toc(hpipe);
        else
            sz = Params.neuron_size;
            Fsi_new = Params.Fsi_new;
            [roi, sig, seeds, bg, bgf, datasmth0, cutoff0, pkcutoff0] = manual_seeds_select(reg, Fsi_new, sz);
        end

        %% parameter init %%
        hpipe = tic;
        [P, options] = par_init(reg);

        %% refine roi %%
        noise = P.sn;
        ispara = Params.refine_roi_ispara;
        [roirf, bgrf, sigupdt, seedsupdt, datasmthf1, cutofff1, pkcutofff1] = refine_roi(reg, sig, bgf, roi, seeds, noise, datasmth0, cutoff0, pkcutoff0, ispara);

        %% refine sig %%
        p = 0; %%% no ar model used %%%
        [~, ~, Puse] = refine_sig(reshape(reg, pixh * pixw, nf), roirf, bgrf, sigupdt, bgf, p, options);

        %% merge roi %%
        corrthres = Params.merge_roi_corrthres;
        [roimrg, sigmrg, seedsmrg, datasmthf2, cutofff2, pkcutofff2] = merge_roi(reg, roirf, sigupdt, seedsupdt, datasmthf1, cutofff1, pkcutofff1, corrthres);

        %% refine roi again %%
        Puse.p = 0;
        Puse.options = options;
        Puse.noise = noise;
        ispara = Params.refine_roi_ispara;
        [roifn, bgfn, sigupdt2, seedsfn] = refine_roi(reg, sigmrg, bgf, roimrg, seedsmrg, Puse.noise, datasmthf2, cutofff2, pkcutofff2, ispara);

        %% refine sig again for raw sig %%
        Puse.p = 0; %%% 0 ar model used %%%
        [sigfnr, ~, ~] = refine_sig(reshape(reg, pixh * pixw, nf), roifn, bgfn, sigupdt2, bgf, Puse.p, Puse.options);
        sigfnr = max(roifn, [], 1)' .* sigfnr;
        roifnr = roifn ./ max(roifn, [], 1);

        %% refine sig again %%
        Puse.p = 2; %%% 2nd ar model used %%%
        [sigfn, bgffn, ~] = refine_sig(reshape(reg, pixh * pixw, nf), roifn, bgfn, sigupdt2, bgf, Puse.p, Puse.options);

        %% save data %%
        imax = max(reg, [], 3);
        file_name_to_save = [path_name, file_base{i}, '_data_processed.mat'];
        file_name_reg = [path_name, file_base{i}, '_reg.mat'];
        
        if ismc
            save(file_name_to_save, 'roifn', 'sigfn', 'seedsfn', 'bgfn', 'bgffn', 'roifnr', 'sigfnr', 'imax', 'pixh', 'pixw', 'corr_score', 'raw_score', 'Params');
        else
            save(file_name_to_save, 'roifn', 'sigfn', 'seedsfn', 'bgfn', 'bgffn', 'roifnr', 'sigfnr', 'imax', 'pixh', 'pixw', 'Params');
        end
        
        if isvis
            save(file_name_to_save, 'imaxn', 'imaxy', '-append');
            save(file_name_reg, 'data', '-v7.3')
        else
            save(file_name_reg, 'reg', '-v7.3')
        end
        
        time2 = toc(hpipe);
        disp(['Done all, total time: ', num2str(time1 + time2), ' seconds'])
    end
end

function min1pipe_init
% parse path, and install cvx if not
%   Jinghao Lu, 11/10/2017

    %%% prepare main folder %%%
    pathname = mfilename('fullpath');
    mns = mfilename;
    lname = length(mns);
    pathtop1 = pathname(1: end - lname);
    
    %%% check if on path %%%
    pathCell = regexp(path, pathsep, 'split');
    if ispc  % Windows is not case-sensitive
        onPath = any(strcmpi(pathtop1(1: end - 1), pathCell)); %%% get rid of filesep %%%
    else
        onPath = any(strcmp(pathtop1(1: end - 1), pathCell));
    end
    
    %%% set path and setup cvx if not on path %%%
    if ~onPath
        pathall = genpath(pathtop1);
        addpath(pathall)
        pathcvx = [pathtop1, 'utilities', filesep, 'cvx', filesep, 'cvx_setup.m'];
        run(pathcvx)
    end
end






