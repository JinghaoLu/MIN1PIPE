function [se, spatialr] = auto_detect_params(a)
    aa = imgaussfilt(a, 2) - imgaussfilt(a, min(size(a)) / 5);
    tb = mean(aa(:));
    aa = max(tb, aa);
    
    %%% focus on the center %%%
    ref = fspecial('gaussian', size(aa), min(size(aa)));
    mask = ref > (max(ref(:)) + min(ref(:))) / 2;
    aat = aa .* mask;
    aat(~mask) = tb;
    aaa = imgaussfilt(aat, 2);
    
    scl = 5;
    b = sort(aaa(:))';
    bb = b - scl * linspace(min(aaa(:)), max(aaa(:)), numel(aaa));
    [~, id] = min(bb);
    thres = b(id);
    bbb = imopen(aaa > thres, strel('disk', 2));
    [l, n] = bwlabeln(bbb);
    
    mx = imregionalmax(aaa);
    idm = find(mx);
    tm = aaa(idm);
    [~, id] = sort(tm, 'descend');
    stp = 5;
    iduse = idm(id(1: stp));
    
    s = zeros(1, stp);
    for i = 1: stp
        ii = l(iduse(i));
        tmp = l == ii;
        st = sum(tmp(:));
        s(i) = sqrt(st / pi);
    end
    
    sf = ceil(mean(s));
    if sf > 6
        se = 5;
        spatialr = 5 / sf;
    else
        spatialr = 1;
        se = sf;
    end
end