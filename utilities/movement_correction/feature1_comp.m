function [mxall, nfrefs] = feature1_comp(m, stt, stp, idbatch)
% compute first feature map of the images
%   Jinghao Lu, 01/26/2018

    [pixh, pixw, nf] = size(m, 'reg');
    mxall = zeros(pixh, pixw, length(stt));
    nfrefs = zeros(1, length(stt));
    nbatch = length(idbatch) - 1;
    for i = 1: nbatch
        Ypara = cell(1, idbatch(i + 1) - idbatch(i));
        for ii = idbatch(i) + 1: idbatch(i + 1)
            Ypara{ii - idbatch(i)} = m.reg(1: pixh, 1: pixw, stt(ii): stp(ii));
        end
        
        stof = idbatch(i);
        parfor ii = idbatch(i) + 1: idbatch(i + 1)
            Ytmp = Ypara{ii - stof};
            [mxall(:, :, ii), nfrefs(ii)] = ref_select(Ytmp);
        end
   end
end





