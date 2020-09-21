function result = rdp(Points, epsilon)
% RDP algorithm for computing key turning points
%   by Jinghao Lu, adapted from Reza Ahmadzadeh

    dmax = 0;
    edx = length(Points);
    for ii = 2: edx - 1
        d = penDistance(Points(:, ii), Points(:, 1), Points(:, edx));
        if d > dmax
            idx = ii;
            dmax = d;
        end
    end
    if dmax > epsilon
        % recursive call
        recResult1 = rdp(Points(:, 1: idx), epsilon);
        recResult2 = rdp(Points(:, idx: edx), epsilon);
        result = [recResult1(:, 1: length(recResult1) - 1), recResult2(:, 1: length(recResult2))];
    else
        result = [Points(:, 1), Points(:, edx)];
    end
    % If max distance is greater than epsilon, recursively simplify
end
    
function d = penDistance(Pp, P1, P2)
    % find the distance between a Point Pp and a line segment between P1, P2.
    d = abs((P2(2, 1) - P1(2, 1)) * Pp(1, 1) - (P2(1, 1) - P1(1, 1)) * Pp(2, 1) + P2(1, 1) * P1(2, 1) - P2(2, 1) * P1(1, 1)) ...
        / sqrt((P2(2, 1) - P1(2, 1)) ^ 2 + (P2(1, 1) - P1(1, 1)) ^ 2);
end
