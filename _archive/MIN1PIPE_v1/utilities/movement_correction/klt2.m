function [d, TPoints, SInliers, isFound] = klt2(S, T, biderr, mq, winsize, oldpoints)
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
    pointTracker = vision.PointTracker('MaxBidirectionalError', biderr, 'NumPyramidLevels', 4, 'BlockSize', winsize);
    initialize(pointTracker, oldpoints, S);
    [points, isFound] = step(pointTracker, T);
    TPoints = points(isFound, :);
    SInliers = oldpoints(isFound, :);
    d = TPoints - SInliers;
end