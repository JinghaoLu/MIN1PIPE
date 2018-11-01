function [inliercur, inlierold, xform] = klt_geo(pold, pcur, hthres, flag)
% get the affine transformation based on point pairs from KLT tracker
%   Jinghao Lu, 05/14/2017

    if nargin < 3 || isempty(hthres)
        hthres = 1.5; %%% matlab default %%%
    end
    
    if nargin < 4
        flag = 1;
    end
        
    try
        if flag == 1
            [xform, inliercur, inlierold] = ...
                estimateGeometricTransform(pcur, pold, ...
                'affine', 'MaxDistance', hthres);
        else
            xform = fitgeotrans(pcur, pold, 'affine');
            inliercur = pcur;
            inlierold = pold;
        end
    catch
        xform = [];
        inliercur = [];
        inlierold = [];
    end
end