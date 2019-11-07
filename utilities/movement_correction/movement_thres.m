function thres = movement_thres(acorr)
    [ct, ctr] = hist(acorr, round(length(acorr) / 10));
    [~, ids] = max(ct);
    y = ct(ids: end);
    y = y(:);
    x = ctr(ids: end);
    x = x(:);
    ys = smooth(y, max(1, round(length(y) / 20)));
    
    yy = linspace(4 * y(1), y(end), length(x));
    tmp = ys(:) - yy(:);
    [~, idx] = min(tmp);
    thres = ctr(ids + idx - 1);
end