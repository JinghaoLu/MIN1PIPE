function mxpsig = feature2_comp(mxall, flag, mag, denom, sz)
% compute 2nd features of the images
%   Jinghao Lu, 01/26/2018

    if nargin < 2 || isempty(flag)
        flag = 0;
    end
    
    if nargin < 3 || isempty(mag)
        mag = 40;
    end
    
    if nargin < 4 || isempty(denom)
        denom = 10;
    end
    
    if nargin < 5 || isempty(sz)
        sz = 9;
    end
    
    mxp = mxall;
    n = size(mxall, 3);
    
    %%% if do anidenoise %%%
    if flag == 1
        parfor i = 1: n
            tmp = anidenoise(imerode(mxp(:, :, i), strel('disk', 1)), sz);
            mxp(:, :, i) = tmp(sz + 1: end - sz, sz + 1: end - sz);
        end
    end
    
    %%% 2nd feature, demons_prep %%%
    mxpsig = mxp; %%% mxpsig: just slightly enhance %%%
    for i = 1: n
        basecur = mxpsig(:, :, i);
        mxpsig(:, :, i) = demons_prep(basecur, mag, denom);
    end
    mxpsig = bsxfun(@minus, mxpsig, min(min(mxpsig, [], 1), [], 2));
end