ds = floor(pixw / 10);
ns = ceil(pixw / ds);
for i = 1: ns
    tmp = m.frame_all(1: pixh, (i - 1) * ds + 1: min(i * ds, pixw), 1: nf);
    tmp1 = reshape(tmp, [], nf);
    tmp2 = diff(tmp1, 1, 2);
    tmp3 = abs(tmp2);
    ids = find(max(tmp1, [], 2) > 0.9 & (max(tmp2, [], 2) + min(tmp2, [], 2)) < 0.05 & max(tmp3, [], 2) > 0.5);
    [idy, idx] = ind2sub([size(tmp, 1), size(tmp, 2)], ids);
    for j = 1: length(idx)
        m.frame_all(idy(j), (i - 1) * ds + idx(j), 1: nf) = 0;
    end
end