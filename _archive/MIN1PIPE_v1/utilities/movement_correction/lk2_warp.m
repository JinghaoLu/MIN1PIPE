function imgreg = lk2_warp(T, PN, p0)
% warp of the LK tracker
%   Jinghao Lu, 05/11/2017

    Ht = eye(3);
    Ht(3, 1: 2) = -fliplr(squeeze(PN - p0)');
    imgreg = imwarp(T, affine2d(Ht), 'OutputView', imref2d(size(T)));
end