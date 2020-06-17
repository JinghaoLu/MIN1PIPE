function [file_name_to_save, filename_raw, filename_reg] = min1pipe_HPC(Fsi, Fsi_new, spatialr, se, ismc, flag, path_name, file_name)
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
    
    if nargin < 7 || isempty(flag)
        [file_name, path_name] = uigetfile('*', 'Select coordinates file', 'MultiSelect', 'on');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%% parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% user defined parameters %%%                                     %%%
    Params.Fsi = Fsi;                                                   %%%
    Params.Fsi_new = Fsi_new;                                           %%%
    Params.spatialr = spatialr;                                         %%%
    Params.neuron_size = se; %%% half neuron size; 9 for Inscopix and 5 %%%
                            %%% for UCLA, with 0.5 spatialr separately  %%%
                                                                        %%%
    %%% fixed parameters (change not recommanded) %%%                   %%%
    Params.anidenoise_iter = 4;                   %%% denoise iteration %%%
    Params.anidenoise_dt = 1/7;                   %%% denoise step size %%%
    Params.anidenoise_kappa = 0.5;       %%% denoise gradient threshold %%%
    Params.anidenoise_opt = 1;                %%% denoise kernel choice %%%
    Params.anidenoise_ispara = 1;             %%% if parallel (denoise) %%%   
    Params.bg_remove_ispara = 1;    %%% if parallel (backgrond removal) %%%
    Params.mc_scl = 0.004;      %%% movement correction threshold scale %%%
    Params.mc_sigma_x = 5;  %%% movement correction spatial uncertainty %%%
    Params.mc_sigma_f = 10;    %%% movement correction fluid reg weight %%%
    Params.mc_sigma_d = 1; %%% movement correction diffusion reg weight %%%
    Params.pix_select_sigthres = 0.8;     %%% seeds select signal level %%%
    Params.pix_select_corrthres = 0.6; %%% merge correlation threshold1 %%%
    Params.refine_roi_ispara = 1;          %%% if parallel (refine roi) %%%
    Params.merge_roi_corrthres = 0.9;  %%% merge correlation threshold2 %%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% get dataset info %%
    [path_name, file_base, file_fmt] = data_info_HPC(file_name, path_name);
    
    hpipe = tic;
    for i = 1: length(file_base)
        
        %%% judge whether do the processing %%%
        filecur = [path_name, file_base{i}, '_data_processed.mat'];
        msg = 'Redo the analysis? (y/n)';
        overwrite_flag = judge_file(filecur, msg);
        
        if overwrite_flag
            %% data cat %%
            Fsi = Params.Fsi;
            Fsi_new = Params.Fsi_new;
            spatialr = Params.spatialr;
            [m, filename_raw, imaxn, imeanf, pixh, pixw, nf] = data_cat(path_name, file_base{i}, file_fmt{i}, Fsi, Fsi_new, spatialr);
            
            %% neural enhancing batch version %%
            filename_reg = [path_name, file_base{i}, '_reg.mat'];
            [m, imaxy, overwrite_flag] = neural_enhance(m, filename_reg, Params);
            
            %% neural enhancing postprocess %%
            if overwrite_flag
                m = noise_suppress(m, imaxy);
            end
            
            %% get rough roi domain %%
            mask = dominant_patch(imaxy);

            %% frame register %%
            if ismc && overwrite_flag
                pixs = min(pixh, pixw);
                Params.mc_pixs = pixs;
                Fsi_new = Params.Fsi_new;
                scl = Params.neuron_size / (7 * pixs);
                sigma_x = Params.mc_sigma_x;
                sigma_f = Params.mc_sigma_f;
                sigma_d = Params.mc_sigma_d;
                se = Params.neuron_size;
                [m, corr_score, raw_score, scl] = frame_reg(m, imaxy, se, Fsi_new, pixs, scl, sigma_x, sigma_f, sigma_d);
                Params.mc_scl = scl; %%% update latest scl %%%
                
                file_name_to_save = [path_name, file_base{i}, '_data_processed.mat'];
                if exist(file_name_to_save, 'file')
                    delete(file_name_to_save)
                end
                save(file_name_to_save, 'corr_score', 'raw_score', '-v7.3');
            end
            
            time1 = toc(hpipe);
            
            %% select pixel %%
            if flag == 1
                sz = Params.neuron_size;
                Fsi_new = Params.Fsi_new;
                sigthres = Params.pix_select_sigthres;
                corrthres = Params.pix_select_corrthres;
                [roi, sig, bg, bgf, seeds, datasmth0, cutoff0, pkcutoff0] = pix_select(m, mask, sz, Fsi_new, sigthres, corrthres);
                time1 = toc(hpipe);
            else
                sz = Params.neuron_size;
                Fsi_new = Params.Fsi_new;
                [roi, sig, seeds, bg, bgf, datasmth0, cutoff0, pkcutoff0] = manual_seeds_select(m, Fsi_new, sz);
            end
            
            %% parameter init %%
            hpipe = tic;
            [P, options] = par_init(m);
            
            %% refine roi %%
            noise = P.sn;
            ispara = Params.refine_roi_ispara;
            [roirf, bgr, sigupdt, seedsupdt, datasmthf1, cutofff1, pkcutofff1] = refine_roi(m, sig, bgf, roi, seeds, noise, datasmth0, cutoff0, pkcutoff0, ispara);
            
            %% refine sig %%
            p = 0; %%% no ar model used %%%
            [sigrf, bgrf, Puse] = refine_sig(m, roirf, bgr, sigupdt, bgf, p, options);
            
            %% merge roi %%
            corrthres = Params.merge_roi_corrthres;
            [roimrg, sigmrg, seedsmrg, datasmthf2, cutofff2, pkcutofff2] = merge_roi(m, roirf, sigrf, seedsupdt, datasmthf1, cutofff1, pkcutofff1, corrthres);
            
            %% refine roi again %%
            Puse.p = 0;
            Puse.options = options;
            Puse.noise = noise;
            ispara = Params.refine_roi_ispara;
            [roifn, bgfn, sigupdt2, seedsfn] = refine_roi(m, sigmrg, bgrf, roimrg, seedsmrg, Puse.noise, datasmthf2, cutofff2, pkcutofff2, ispara);
            
            %% refine sig again for raw sig %%
            Puse.p = 0; %%% 0 ar model used %%%
            [sigfnr, ~, ~] = refine_sig(m, roifn, bgfn, sigupdt2, bgf, Puse.p, Puse.options);
            sigfnr = max(roifn, [], 1)' .* sigfnr;
            roifnr = roifn ./ max(roifn, [], 1);
            
            %% refine sig again %%
            Puse.p = 2; %%% 2nd ar model used %%%
            Puse.options.p = 2;
            [sigfn, bgffn, ~, spkfn] = refine_sig(m, roifn, bgfn, sigupdt2, bgf, Puse.p, Puse.options);
            sigfn = max(roifn, [], 1)' .* sigfn;
            roifn = roifn ./ max(roifn, [], 1);
%             dff = sigfn ./ mean(sigfn, 2);
            dff = sigfn ./ mean(bgffn);
            
            %% save data %%
            stype = parse_type(class(m.reg(1, 1, 1)));
            nsize = pixh * pixw * nf * stype; %%% size of single %%%
            nbatch = batch_compute(nsize);
            ebatch = ceil(nf / nbatch);
            idbatch = [1: ebatch: nf, nf + 1];
            nbatch = length(idbatch) - 1;
            imax = zeros(pixh, pixw);
            for j = 1: nbatch
                tmp = m.reg(1: pixh, 1: pixw, idbatch(j): idbatch(j + 1) - 1);
                imax = max(cat(3, max(tmp, [], 3), imax), [], 3);
            end
            
            file_name_to_save = [path_name, file_base{i}, '_data_processed.mat'];
            if exist(file_name_to_save, 'file')
                if ismc
                    load(file_name_to_save, 'raw_score', 'corr_score')
                end
                delete(file_name_to_save)
            end
            
            if ismc
                save(file_name_to_save, 'roifn', 'sigfn', 'dff', 'seedsfn', 'spkfn', 'bgfn', 'bgffn', 'roifnr', 'sigfnr', 'imax', 'pixh', 'pixw', 'corr_score', 'raw_score', 'Params', '-v7.3');
            else
                save(file_name_to_save, 'roifn', 'sigfn', 'dff', 'seedsfn', 'spkfn', 'bgfn', 'bgffn', 'roifnr', 'sigfnr', 'imax', 'pixh', 'pixw', 'Params', '-v7.3');
            end
            
            save(file_name_to_save, 'imaxn', 'imaxy', '-append');
            time2 = toc(hpipe);
            disp(['Done all, total time: ', num2str(time1 + time2), ' seconds'])
        else
            filename_raw = [path_name, file_base{i}, '_frame_all.mat'];
            filename_reg = [path_name, file_base{i}, '_reg.mat'];
            file_name_to_save = filecur;
            
            time2 = toc(hpipe);
            disp(['Done all, total time: ', num2str(time2), ' seconds'])
        end
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
        cvx_dir = [pathtop1, 'utilities'];
        if ~exist([cvx_dir, filesep, 'cvx'], 'dir')
            if ispc
                cvxl = 'http://web.cvxr.com/cvx/cvx-w64.zip';
            elseif isunix
                cvxl = 'http://web.cvxr.com/cvx/cvx-a64.zip';
            elseif ismac
                cvxl = 'http://web.cvxr.com/cvx/cvx-maci64.zip';
            end
            disp('Downloading CVX');
            unzip(cvxl, cvx_dir);
        end
        pathcvx = [cvx_dir, filesep, 'cvx', filesep, 'cvx_setup.m'];
        run(pathcvx)
    end
end






