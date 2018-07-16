function savef(filename, flag, varargin)
% save to disk with appending capability using uncompressed h5write
%   Jinghao Lu 06/21/2018

    %%% Extract the variable values %%%
    vars = cell(size(varargin));
    for i = 1:numel(vars)
        vars{i} = evalin('caller', varargin{i});
    end

    %%% Separate numeric arrays from the rest %%%
    isnum = cellfun(@(x) isa(x, 'numeric'), vars);

    %%% Append .mat if necessary %%%
    [filepath, filebase, ext] = fileparts(filename);
    if isempty(ext)
        filename = fullfile(filepath, [filebase '.mat']);
    end

    if all(isnum)
        %%% Save a dummy variable, just to create the file %%%
        if ~exist(filename, 'file')
            dummy = 0;
            save(filename, '-v7.3', 'dummy');

            %%% Delete the dummy, if necessary, just in case the user supplied a
            %%% variable called dummy %%%
            fid = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
            H5L.delete(fid,'dummy','H5P_DEFAULT');
            H5F.close(fid);
        end
    else
        s = struct;
        for i = 1:numel(isnum)
            if ~isnum(i)
                s.(varargin{i}) = vars{i};
            end
        end
        save(filename, '-v7.3', '-struct', 's');
    end


    %%% Save all numeric variables %%%
    if flag == 1 %%% save independent variables %%%
        for i = 1:numel(isnum)
            if ~isnum(i)
                continue
            end
            varname = ['/', varargin{i}];
            h5create(filename, varname, size(vars{i}), 'DataType', class(vars{i}));
            h5write(filename, varname, vars{i});
        end
    else %%% save concatenated variables %%%
%         h5create(filename1, varname, [20, 20, Inf], 'ChunkSize',[5, 5, 5], 'DataType', 'single');
        m = matfile(filename);
        for i = 1:numel(isnum)
            if ~isnum(i)
                continue
            end
            varname = ['/', varargin{i}];
            tmp = whos(m);
            ids = [];
            for ii = 1: length(tmp)
                if strcmp(varargin{i}, tmp(ii).name)
                    ids = ii;
                    break
                end
            end
            
            if ~isempty(ids)
                sz = tmp(ids).size;
                stt = sz(end) + 1;
            else
                sz = [size(vars{i}), 1];
                sz = sz(1: 3); %%% force to concatenate on 3rd dimension %%%
                stt = 1;
            end
            
            try
                h5create(filename, varname, [sz(1: end - 1), Inf], 'ChunkSize', [sz(1: end - 1), 1], 'DataType', class(vars{i}));
            catch
                
            end
            start = [ones(1, length(sz) - 1), stt];
            count = [size(vars{i}), 1];
            count = count(1: 3);
            h5write(filename, varname, vars{i}, start, count);
        end
    end
end
