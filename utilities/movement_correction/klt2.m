function [d, TPoints, SInliers, isFound] = klt2(S, T, biderr, mq, winsize, oldpoints, maskc)
% KLT tracker between the two images
%   Jinghao Lu, 08/25/2017

    if nargin < 3 || isempty(biderr)
        biderr = 2;
    end
    if nargin < 4 || isempty(mq)
        mq = 0.01;
    end
    if nargin < 5 || isempty(winsize)
        winsize = [31, 31];
    end
    if nargin < 6 || isempty(oldpoints)
        oldpoints = detectMinEigenFeatures(S, 'MinQuality', mq);
        oldpoints = oldpoints.Location;        
    end
    flagmask = true;
    if nargin < 7 || isempty(maskc)
        flagmask = false;
    end
    pointTracker = vision.PointTracker('MaxBidirectionalError', biderr, 'NumPyramidLevels', 4, 'BlockSize', winsize);
    initialize(pointTracker, oldpoints, S);
    [points, isFound] = step(pointTracker, T);
    TPoints = points(isFound, :);
    SInliers = oldpoints(isFound, :);
    if flagmask
        idS = [];
        for i = 1: size(SInliers, 1)
            if maskc(round(SInliers(i, 2)), round(SInliers(i, 1)))
                idS = [idS, i];
            end
        end
        TPoints = TPoints(idS, :);
        SInliers = SInliers(idS, :);
    end
    d = TPoints - SInliers;
end