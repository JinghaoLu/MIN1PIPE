function [scrout, imgout, xformout] = klt_ref_track(Y, imref, nmaxloop, maskc)
% track with KLT tracker relative to imref
%   Jinghao Lu, 02/20/2018
    
    if nargin < 3 || isempty(nmaxloop)
        nmaxloop = 10;
    end
    
    if nargin < 4 || isempty(maskc)
        maskc = true(size(imref));
    end
        
    %% KLT tracker %%
    %%% simple version: no loop just once %%%
    nf = size(Y, 3);
    scrout = zeros(1, nf);
    xformout = cell(1, nf);
    imgout = zeros(size(Y));
    
    %%% loop through all neighboring frame pairs %%%
    for i = 1: nf
        imcur = Y(:, :, i); 
        
        xforms = cell(size(imref, 3), nmaxloop + 1);
        imgs = cell(size(imref, 3), nmaxloop + 1);
        scrs = zeros(size(imref, 3), nmaxloop + 2); %%% additional score for untransformed %%%
            
        %%% nmaxloop times realizations %%%
        for ii = 1: nmaxloop
            [img, xf] = klt2_reg(imref, imcur, 1, maskc);
            scr = get_trans_score_ref(img, imref, maskc);
            xforms{ii} = xf;
            imgs{ii} = img;
            scrs(ii) = scr;
        end
        
        %%% gather additional fitgeotrans and scores of unregistered frame pairs %%%
        [img, xf] = klt2_reg(imref, imcur, 2, maskc);
        scr = get_trans_score_ref(img, imref, maskc);
        xforms{nmaxloop + 1} = xf;
        imgs{nmaxloop + 1} = img;
        scrs(nmaxloop + 1) = scr;
        scr = get_trans_score_ref(img, imref, maskc);
        scrs(nmaxloop + 2) = scr;
        
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
        imgout(:, :, i) = img;
        xformout{i} = xform;
    end
end