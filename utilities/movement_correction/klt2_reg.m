function [img, xform] = klt2_reg(imref, imcur, flag, maskc)
% track and register 2 images with KLT tracker
%   Jinghao Lu, 11/15/2017

    if nargin < 3
        flag = 1;
    end
    
    if nargin < 4 || isempty(maskc)
        maskc = true(size(imref));
    end
    
    [~, pcur, pold, ~] = klt2(imref, imcur, [], [], [], [], maskc);
    [~, ~, xform] = klt_geo(pold, pcur, [], flag); %%% 1: regular RANSAC trans; else: fitgeotrans %%%
    if ~isempty(xform)
        img = klt_warp(imcur, xform);
    else
        img = imcur;
        xform = affine2d(diag(ones(1, 3)));
    end
end
