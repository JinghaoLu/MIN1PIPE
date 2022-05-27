function [m_out, imaxn] = remove_dp(m_in, vname)
    [pixh, pixw, nf] = size(m_in, vname);
    m_in.Properties.Writable = true;
%     [X, Y] = meshgrid(1: pixw, 1: pixh);
%     X = X(:);
%     Y = Y(:);
%     ref = [X, Y];
%     ctr = round([pixw, pixh] / 2);
%     dis = vecnorm(ref - ctr, 2, 2);
%     dis = reshape(dis, pixh, pixw);
%     dis = max(dis(:)) - dis;
%     kn = 1 ./ (1 + exp(-0.5 * (dis - 0.01 * max(ctr))));
    kn = normalize(fspecial('gaussian', [pixh, pixw], max(pixh, pixw)));
    kn1 = kn < 0.02;
    kn2 = kn < 0.03;
    idt = xor(kn1, kn2);
    [l, n] = bwlabeln(idt);
    [l1, n1] = bwlabeln(kn1);
    
    eval(['stype = class(m_in.', vname, '(1, 1, 1));'])
    stype = parse_type(stype);
    nsize = pixh * pixw * nf * stype; %%% size of single %%%
    nbatch = batch_compute(nsize);
    ebatch = ceil(nf / nbatch);
    
    idbatch = [1: ebatch: nf, nf + 1];
    imx1 = 0;
    imn1 = inf;
    imaxn = zeros(pixh, pixw);
    
    for ib = 1: nbatch
        eval(['frame_all = m_in.', vname, '(1: pixh, 1: pixw, idbatch(ib): idbatch(ib + 1) - 1);'])
        frame_all = reshape(frame_all, pixh * pixw, []);
        for j = 1: n
            idtt = l == j;
            tmp2 = frame_all(idtt(:), :);
            mtmp = median(tmp2, 1);  
            idtt = l1 == j;
            frame_all(idtt(:), :) = repmat(mtmp, sum(idtt(:)), 1);
        end
        frame_all = reshape(frame_all, pixh, pixw, []);
        imx1 = max(imx1, max(max(max(frame_all))));
        imn1 = min(imn1, min(min(min(frame_all))));
        imaxn = max(cat(3, imaxn, max(frame_all, [], 3)), [], 3);
        eval(['m_in.', vname, '(1: pixh, 1: pixw, idbatch(ib): idbatch(ib + 1) - 1) = frame_all;'])
    end
    
    m_out = normalize_batch(m_in.Properties.Source, 'frame_allt', imx1, imn1, idbatch);
end