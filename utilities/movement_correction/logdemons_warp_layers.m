function [mxout] = logdemons_warp_layers(mxin, xffn, ldfn)
% warp hierarchical layers with affine transforms and motion field
%   Jinghao Lu, 01/02/2018

    mxtmp = mxin;
    n = length(xffn);
    parfor i = 1: n
        imgs = mxtmp{i};
        nfc = size(imgs, 3);
        if ~isempty(xffn{i})
            ntemp = length(xffn{i});
            for j = 1: ntemp
                for iframe = 1: nfc
                    imgs(:, :, iframe) = klt_warp(imgs(:, :, iframe), xffn{i}{j});
                end
                nldtmp = length(ldfn{i}{j});
                for k = 1: nldtmp
                    ldtemp = ldfn{i}{j}{k};
                    for iframe = 1: nfc
                        imgs(:, :, iframe) = iminterpolate(imgs(:, :, iframe), ldtemp(:, :, 1), ldtemp(:, :, 2));
                    end
                end
            end
        end
        mxtmp{i} = imgs;
    end
    mxout = mxtmp;
end