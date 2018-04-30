function fps = feature_select(mxp)
% find salient regional max as feature points
%   Jinghao Lu, 08/25/2017

    fps = cell(1, size(mxp, 3));
    for ii = 1: size(mxp, 3)
        tmp = mxp(:, :, ii);
        rmax = imregionalmax(tmp);
%         ithres = prctile(tmp(:), 98);
        ithres = (max(tmp(:)) + min(tmp(:))) / 3;
        rmt = rmax .* tmp > ithres;
        [idh, idw] = find(rmt);
        fps(ii) = {[idh, idw]};
        [~, idt] = sort(tmp(rmt > 0), 'descend');
        fps{ii} = fps{ii}(idt, :);
    end
end
