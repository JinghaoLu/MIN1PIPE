function s = auto_detect_params(a)
    aa = a - imgaussfilt(a, min(size(a)) / 5);
    aa = max(mean(aa(:)), aa);
    aaa = imgaussfilt(aa,2);
    scl = 5;
    b = sort(aaa(:))';
    bb = b - scl * linspace(min(aaa(:)), max(aaa(:)), numel(aaa));
    [~, id] = min(bb);
    thres = b(id);
    bbb = imopen(aaa > thres, strel('disk', 2));
    [l, n] = bwlabeln(bbb);
    s = zeros(1, n);
    for i = 1: n
        tmp = l == i;
        st = sum(tmp(:));
        s(i) = sqrt(st / pi);
    end
end