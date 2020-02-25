function [PN, imgo, scr] = lk_ref_track(Y, imref, maskc)
% track image cluster Y with LK tracker relative to the reference frame
%   Jinghao Lu, 01/30/2018

    if nargin < 3 || isempty(maskc)
        maskc = true(size(imref));
    end

    %%% track frames %%%
    n = size(Y, 3);
    PN = zeros(2, n);
    isf = zeros(1, n);
    scr = zeros(1, n);
    imgo = Y;
    for i = 1: n
        imcur = Y(:, :, i);
        scr(i) = get_trans_score_ref(imcur, imref, maskc);
        [p0, winsize] = cen_gen(cat(3, imref, imcur));
        [temp, isf(i), ~] = lk2_track(imref, imcur, p0, winsize);
        PN(:, i) = temp - p0;
        if ~isf(i)
            imreft = feature2_comp(imref);
            imcurt = feature2_comp(imcur);
            [p0, winsize] = cen_gen(cat(3, imreft, imcurt));
            [temp, isf(i), ~] = lk2_track(imreft, imcurt, p0, winsize);
            PN(:, i) = temp - p0;
            if ~isf(i)
                [temp, ~] = imregdft(imref, imcur);
                PN(:, i) = temp;
            end
        end
        p0 = [0; 0];
        imgo(:, :, i) = lk2_warp(imcur, PN(:, i), p0);
        
        if sqrt(PN(1, i) ^ 2 + PN(2, i) ^ 2) > 4 * scr(i)
            PN(:, i) = p0;
            imgo(:, :, i) = imcur;
        end
    end
end