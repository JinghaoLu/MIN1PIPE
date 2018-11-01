function mask = roi_domain(frame)
% get dominant roi domain
%   Jinghao Lu 01/16/2016

%     mx = mean(frame, 3);
%     mask = mx > prctile(mx(:), 25);
    mx = mean(frame, 3);
    mask = normalize(imgaussfilt(mx, max(size(mx)) / 5)) > 0.25;
%     mask = normalize(imgaussfilt(feature2_comp(mx, 0, 200, 5), size(mx) / 3)) > 0.5;
    [l, n] = bwlabeln(mask);
    s = zeros(1, n);
    for i = 1: n
        s(i) = sum(sum(l == i));
    end
    [~, idt] = max(s);
    mask = l == idt;
end