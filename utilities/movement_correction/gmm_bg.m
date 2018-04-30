function thres = gmm_bg(acorr)
% get background threshold
%   Jinghao Lu, 05/16/2016

    G = fitgmdist(acorr(:), 2);
    [~, iuse] = min(G.mu);
    scall = posterior(G, acorr(:));
    idx = scall(:, iuse) > 0.5;
    tmp = acorr(idx);
    thres = prctile(tmp, 99);
end