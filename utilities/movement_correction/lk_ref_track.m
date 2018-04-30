function [PN, imgo, scr] = lk_ref_track(Y, imref)
% track image cluster Y with LK tracker relative to the reference frame
%   Jinghao Lu, 01/30/2018

    %%% track frames %%%
    n = size(Y, 3);
    PN = zeros(2, n);
    isf = zeros(1, n);
    scr = zeros(1, n);
    imgo = Y;
    for i = 1: n
        imcur = Y(:, :, i);
        [p0, winsize] = cen_gen(cat(3, imref, imcur));
        [temp, isf(i), ~] = lk2_track(imref, imcur, p0, winsize);
        PN(:, i) = temp - p0;
        imgo(:, :, i) = lk2_warp(imcur, temp, p0);
        scr(i) = get_trans_score_ref(imcur, imref);
    end
end