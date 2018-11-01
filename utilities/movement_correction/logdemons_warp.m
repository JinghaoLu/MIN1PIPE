function logdemons_warp(m, sxt, syt, sttnni, stpnni)
    [pixh, pixw, ~] = size(m, 'reg');
    for j = sttnni: stpnni
        tmp = m.reg(1: pixh, 1: pixw, j);
        for jj = 1: length(sxt)
            for k = 1: length(sxt{jj})
                tmp = iminterpolate(tmp, sxt{jj}{k}, syt{jj}{k});
            end
        end
        m.reg(1: pixh, 1: pixw, j) = tmp;
    end
end