function mask = dominant_patch(Y, thres)
% compute the dominant patch of all images, and find the mask
%   Jinghao Lu, 11/10/2017

    if nargin < 2
        thres = 0.1;
    end
    [pixh, pixw, ~] = size(Y);
    t = max(Y, [], 3);
    mskt = double(normalize(demons_prep(t)) > thres);
    msktt = medfilt2(mskt, [ceil(pixh / 50), ceil(pixw / 50)]);
    if ~sum(msktt(:)) == 0
        mskt = msktt;
    end
    bb = regionprops(mskt, 'BoundingBox');
    bb = bb.BoundingBox;
    brg = [ceil(bb(2)), min(pixh, ceil(bb(2)) + bb(4) - 1); ceil(bb(1)), min(pixw, ceil(bb(1)) + bb(3) - 1)];
    ci = regionprops(mskt, 'ConvexImage');
    ci = ci.ConvexImage;
    ciuse = ci(1: diff(brg(1, :)) + 1, 1: diff(brg(2, :)) + 1);
    mskt(brg(1, 1): brg(1, 2), brg(2, 1): brg(2, 2)) = ciuse;
    mskt = normalize(imgaussfilt(mskt, max(pixh, pixw) / 5)) > thres;
    [ll, n] = bwlabeln(mskt);
    ss = zeros(1, n);
    for i = 1: n
        ss(i) = sum(sum(ll == i));
    end
    [~, idmsk] = max(ss);
    mask = ll == idmsk;
end