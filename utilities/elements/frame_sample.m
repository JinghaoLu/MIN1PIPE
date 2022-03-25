function [se, spatialr] = frame_sample(m)
    if contains(m, 'tif') || contains(m, 'tiff')
        f = imfinfo(m);
        n = length(f);
    elseif contains(m, 'avi')
        f = VideoReader;
        n = floor(f.FrameRate * f.Duration);
    elseif contains(m, 'mat')
        f = matfile(m);
        [pixh, pixw, n] = size(f.frame_all);
    end
    nu = min(100, n);
    s = cell(1, nu);
    id = randsample(n, nu);
    for i = 1: nu
        if contains(m, 'tif') || contains(m, 'tiff')
            tmp = double(imread(m, id(i)));
        elseif contains(m, 'avi')
            tmp = double(read(f, id(i)));
        elseif contains(m, 'mat')
            tmp = double(m.frame_all(1: pixh, 1: pixw, id(i)));
        end
        st = auto_detect_params(tmp);
        s{i} = st;
    end
    s = cell2mat(s);
    sf = floor(mean(s));
    if sf > 6
        spatialr = 3 / sf;
        se = 5;
    else
        spatialr = 1;
        se = sf;
    end
end