function img_warp = klt_warp(img, xform)
% warp affine transformation from KLT tracker
%   Jinghao Lu, 05/14/2017

    img_warp = imwarp(img, xform, 'OutputView', imref2d(size(img)));
end