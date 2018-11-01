function mask = roi_domain(frame)
% get dominant roi domain
%   Jinghao Lu 01/16/2016

%     mx = mean(frame, 3);
%     mask = mx > prctile(mx(:), 25);
    mx = mean(frame, 3);
    mask = normalize(imgaussfilt(mx, max(size(mx)) / 5)) > 0.25;
end