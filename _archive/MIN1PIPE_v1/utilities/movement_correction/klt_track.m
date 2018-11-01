function [scrout, imgout, xformout] = klt_track(mxall, flag)
% track with KLT tracker for image cluster, with n-1 xformout
%   Jinghao Lu, 01/30/2018

    if nargin < 2
        flag = 1;
    end
    
    %% KLT tracker %%
    %%% simple version: no loop just once %%%
    nf = size(mxall, 3);
    nmaxloop = 10;
    scrout = zeros(1, nf - 1);
    imgout = cat(3, mxall(:, :, 1), zeros([size(mxall(:, :, 1)), nf - 1]));
    xformout = [];
    
    %%% loop through all neighboring frame pairs %%%
    for i = 1: nf - 1
        imcur = normalize(mxall(:, :, i + 1)); 
%         imcur = im_bg_suppress(imcur);
        imref = normalize(mxall(:, :, i)); %%% raw ref %%%
%         imref = im_bg_suppress(imref);
        
        xforms = cell(size(imref, 3), nmaxloop + 1);
        imgs = cell(size(imref, 3), nmaxloop + 1);
        scrs = zeros(size(imref, 3), nmaxloop + 2); %%% additional score for untransformed %%%
        for k = 1: size(imref, 3)
            imreft = imref(:, :, k);
            
            %%% nmaxloop times realizations %%%
            parfor ii = 1: nmaxloop
                [img, xf] = klt2_reg(imreft, imcur, 1);
                scr = get_trans_score(cat(3, imreft, img), flag);
                xforms{k, ii} = xf;
                imgs{k, ii} = img;
                scrs(k, ii) = scr;
            end
            
            %%% gather additional fitgeotrans and scores of unregistered frame pairs %%%
            [img, xf] = klt2_reg(imreft, imcur, 2);
            scr = get_trans_score(cat(3, imreft, img), flag);
            xforms{k, nmaxloop + 1} = xf;
            imgs{k, nmaxloop + 1} = img;
            scrs(k, nmaxloop + 1) = scr;
            scr = get_trans_score(cat(3, imreft, imcur), flag);
            scrs(k, nmaxloop + 2) = scr;
        end
        
        %%% find the transform to use %%%
        scr = min(scrs(:));
        [x, y] = find(scrs == scr);
        if y(1) < nmaxloop + 2
            img = imgs{x(1), y(1)};
            xform = xforms{x(1), y(1)};
        else
            img = imcur;
            xform = affine2d(diag(ones(1, 3)));
        end
        
        scrout(i) = scr;
        imgout(:, :, i + 1) = img;
        xformout = [xformout, xform];
    end
end