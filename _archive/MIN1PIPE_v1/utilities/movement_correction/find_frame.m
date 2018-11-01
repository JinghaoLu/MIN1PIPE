function [imref, imcur] = find_frame(mxin, idclustfn, idin)
% find original level frame based on the current level index
%   Jinghao Lu, 10/11/2017

    n = length(idclustfn);
    %%%%% imref id %%%%%
    itmp = idin;
    for jj = n: -1: 1
        itmp = idclustfn{jj}{itmp}(end);
    end
    idr = itmp;

    %%%%% imcur id %%%%%
    itmp = idin + 1;
    for jj = n: -1: 1
        itmp = idclustfn{jj}{itmp}(1);
    end
    idc = itmp;

    imref = mxin(:, :, idr);
    imcur = mxin(:, :, idc);
end