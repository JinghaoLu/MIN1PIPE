function [se, spatialr] = frame_sample(m)
    f = imfinfo(m);
    n = length(f);
    nu = min(100, n);
    s = cell(1, nu);
    id = randsample(n, nu);
    for i = 1: nu
        tmp = double(imread(m, id(i)));
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