function [PN, isFound, tt] = lk2_track(S, T, p0, winsize)
% Lucas-Kanade tracker with Matlab built-in functions
%   Jinghao Lu, 05/11/2017

    pointTracker = vision.PointTracker('MaxBidirectionalError', 2, 'BlockSize', winsize, 'NumPyramidLevels', 4);
    initialize(pointTracker, flipud(p0)', S);
    [points, isFound, tt] = step(pointTracker, T);
    PN = fliplr(points)';
end
