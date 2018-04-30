function plot_contour(roifn, sigfn, seedsfn, imax, pixh, pixw)
% plot contour of reconstructed ROIs on imax
%   Jinghao Lu, 06/16/2016

    %%% prepare contours %%%
    [x, y] = ind2sub([pixh, pixw], seedsfn);
    ids = 1: length(x);
    cntrs = cell(1, length(ids));
    thres = 0.8;
    for i = 1: length(ids)
        tmp = full(reshape(roifn(:, ids(i)), pixh, pixw) * max(sigfn(ids(i), :), [], 2));
        tmp = imgaussfilt(tmp, 3);
        lvl = max(max(tmp)) * thres;
        cntrs{i} = contour(flipud(tmp), [lvl, lvl]);
        cntrs{i} = [cntrs{i}(:, 2: end - 1), cntrs{i}(:, 2)];
    end
    
    %%% plot %%%
    tmpp = imax;
    imagesc(tmpp, [0, 0.8])
    hold on
    for i = 1: length(ids)
        plot(cntrs{i}(1, :), pixh - cntrs{i}(2, :), 'r')
        text(y(i), x(i), num2str(i))
    end
    hold off
    title('Neural Contours')
end
