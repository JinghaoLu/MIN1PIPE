function [frame_all, pixh, pixw, nf] = data_cat(path_name, file_base, file_fmt, Fsi, Fsi_new, ratio, ispara)
% Concatinate data pieces from raw 4GB-tiff chunks
%   or avi files from UCLA miniscope
%   parallel or serial version supported
%   Jinghao Lu 01/16/2016

    hcat = tic;
    %% initialization %%
    %%% initialize parameters %%%
    if nargin < 3 || isempty(Fsi)
        defpar = default_parameters;
        Fsi = defpar.Fsi;
    end

    if nargin < 4 || isempty(Fsi_new)
        defpar = default_parameters;
        Fsi_new = defpar.Fsi_new;
    end

    if nargin < 5 || isempty(ratio)
        defpar = default_parameters;
        ratio = defpar.spatialr;
    end
    
    %%% prepare parallel computing %%%
    if ispara
        if isempty(gcp('nocreate'))
            parpool(feature('numCores'));
        end
    end
    
    %%% initialize to get info %%%
    if contains(file_fmt, 'avi')
        dirst = [dir([path_name, file_base, '.avi']); dir([path_name, file_base, '*', '.avi'])];
    elseif contains(file_fmt, 'tiff')
        dirst = [dir([path_name, file_base, '.tiff']); dir([path_name, file_base, '-*', '.tiff'])];
    elseif contains(file_fmt, 'tif')
        dirst = [dir([path_name, file_base, '.tif']); dir([path_name, file_base, '-*', '.tif'])];
    end
    dirs = cell(1, length(dirst));
    nft = zeros(1, length(dirst));
    idx = zeros(1, length(dirst));

    %% get info %%
    disp('Begin collecting datasets info')
    for i = 1: length(dirst)
        dirs{i} = dirst(i).name;
        if contains(file_fmt, 'avi')
            info = VideoReader([path_name, dirs{i}]);
            nft(i) = floor(info.Duration * info.FrameRate);
            temp1 = strfind(dirs{i}, file_base);
            temp2 = strfind(dirs{i}, '.');
            idx(i) = str2double(dirst(i).name(temp1 + length(file_base): temp2 - 1));
        else
            info = imfinfo([path_name, dirs{i}]);
            nft(i) = numel(info);
            temp = strfind(dirst(i).name, '-');
            if isempty(temp)
                idx(i) = 0;
            else
                idx(i) = str2double(dirst(i).name(temp + 1: temp + 3));
            end
        end
    end
    time = toc(hcat);
    disp(['Done collecting datasets info, time: ', num2str(time)])

    %% initialization to get frames %%
    if contains(file_fmt, 'avi')
        info = VideoReader([path_name, dirs{1}]);
    else
        info = imfinfo([path_name, dirs{1}]);
    end
    pixw = info(1,1).Width;
    pixh = info(1,1).Height;
    pixw = round(pixw * ratio);
    pixh = round(pixh * ratio);
    ds = Fsi / Fsi_new;
    [~, induse] = sort(idx);
    dirs = dirs(induse);
    nft = nft(induse);
    nf = sum(nft);
    f_idx = zeros(1, nf);
    dir_idx = cell(1, nf);
    for i = 1: length(nft)
        f_idx(sum(nft(1: i - 1)) + 1: sum(nft(1: i))) = 1: nft(i);
        dir_idx(sum(nft(1: i - 1)) + 1: sum(nft(1: i))) = repmat(dirs(i), 1, nft(i));
    end
    f_use = 1: ds: nf;
    f_idx = f_idx(f_use);
    dir_idx = dir_idx(f_use);
    frame_all = zeros(pixh, pixw, ceil(nf / ds));
    nf = size(frame_all, 3);
    
    %% collect data %%
    disp('Begin data cat')
    if contains(file_fmt, 'avi')
        ft_idx = (f_idx - 1) / Fsi; %%% time stamp of the used frames %%%
        if ispara
            parfor i_frame = 1: length(f_use)
                infot = VideoReader([path_name, dir_idx{i_frame}]);
                infot.CurrentTime = ft_idx(i_frame);
                frame = double(readFrame(infot));
                frame = imresize(frame, [pixh, pixw]);
                frame_all(:, :, i_frame) = frame;
                if mod(i_frame, 100) == 0
                    disp([num2str(i_frame), '/', num2str(nf)])
                end
            end
        else
            dir_old = '';
            for i_frame = 1: length(f_use)
                if ~strcmp(dir_idx{i_frame}, dir_old)
                    infot = VideoReader([path_name, dir_idx{i_frame}]);
                    dir_old = dir_idx{i_frame};
                end
                infot.CurrentTime = ft_idx(i_frame);
                frame = double(readFrame(infot));
                frame = frame(:, :, 1);
                frame = imresize(frame, [pixh, pixw]);
                frame_all(:, :, i_frame) = frame;
                if mod(i_frame, 100) == 0
                    disp([num2str(i_frame), '/', num2str(nf)])
                end
            end
        end
    else
        if ispara
            parfor i_frame = 1: length(f_use)
                frame = double(imread([path_name, dir_idx{i_frame}], f_idx(i_frame)));
                frame = imresize(frame, [pixh, pixw]);
                frame_all(:, :, i_frame) = frame;
                if mod(i_frame, 100) == 0
                    disp([num2str(i_frame), '/', num2str(nf)])
                end
            end
        else
            for i_frame = 1: length(f_use)
                frame = double(imread([path_name, dir_idx{i_frame}], f_idx(i_frame)));
                frame = imresize(frame, [pixh, pixw]);
                frame_all(:, :, i_frame) = frame;
                if mod(i_frame, 100) == 0
                    disp([num2str(i_frame), '/', num2str(nf)])
                end
            end
        end
    end
    time = toc(hcat);
    disp(['Done data cat, time: ', num2str(time)])
end



