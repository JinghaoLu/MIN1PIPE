function imout = demons_prep(img, mag, denom)
% create feature maps good for demons
%   Jinghao Lu, 05/11/2017

    if nargin < 2
        mag = 20;
    end
    if nargin < 3
        denom = 3;
    end
    ithres = (max(img(:)) + min(img(:))) / denom;
    tmp = sigmf(img, [mag, ithres]);
%     tmp = sigmf(img, [mag, prctile(img(:), 98)]);
    tmp = tmp .* (tmp >= median(tmp(:)));
    tmp(tmp == 0) = min(tmp(tmp > 0));
    imout = normalize(tmp);
end