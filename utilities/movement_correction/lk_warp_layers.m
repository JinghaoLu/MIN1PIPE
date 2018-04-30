function [mxout] = lk_warp_layers(mxin, xffn)
% warp each hierarchical layers based on affine transform
%   Jinghao Lu, 09/01/2017

    mxtmp = mxin;
    n = length(xffn);
    for i = 1: n
        imgs = mxtmp{i};
        nfc = size(imgs, 3);
        if ~isempty(xffn{i})
            ntemp = length(xffn{i});
            for j = 1: ntemp
                for iframe = 1: nfc
                    imgs(:, :, iframe) = klt_warp(imgs(:, :, iframe), xffn{i}{j});
                end
            end
        end
        mxtmp{i} = imgs;
    end
    mxout = mxtmp;
end