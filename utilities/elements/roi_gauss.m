function roig = roi_gauss(roi)
% gauss filt for the post process
%  Jinghao Lu 05/20/2017

    [d1, d2, d3] = size(roi);
    d = d1 * d2;
    roig = zeros(d, d3);
    swin = 1;
    cthres = 0.2;
    if ndims(roi) ~= 3
        d3 = size(roi, ndims(roi));
        roi = reshape(roi, d1, d2, d3);
    end
    for i = 1: d3
        roig(:, i) = reshape(normalize(imgaussfilt(full(roi(:, :, i)), swin)) > cthres, d, 1);
    end
end


