function [m_out, Params, pixh, pixw] = downsamp(path_name, file_base, m_in, Params, aflag, imaxn)
    filename = [path_name, file_base, '_frame_all.mat'];
    msg = 'Overwrite raw .mat file (data)? (y/n)';
    overwrite_flag = judge_file(filename, msg);
    if overwrite_flag
        if exist(filename, 'file')
            delete(filename);
        end
        if aflag
            [se, spatialr] = auto_detect_params(imaxn);
            Params.neuron_size = se;
            Params.spatialr = spatialr;
            [m_out, pixh, pixw] = downsamp_unit(m_in, spatialr);
        else
            [m_out, pixh, pixw] = downsamp_unit(m_in, Params.spatialr);
        end
    else
        m_out = matfile(filename, 'Writable', true);
    end
end