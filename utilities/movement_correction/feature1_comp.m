function [mxall, nfrefs] = feature1_comp(Y, stt, stp)
% compute first feature map of the images
%   Jinghao Lu, 01/26/2018

    [pixh, pixw, ~] = size(Y);
    mxall = zeros(pixh, pixw, length(stt));
    nfrefs = zeros(1, length(stt));
    Ypara = cell(1, length(stt));
    for i = 1: length(stt)
        Ypara{i} = Y(:, :, stt(i): stp(i));
    end
    parfor ii = 1: length(stt)
        Ytmp = Ypara{ii};
        [mxall(:, :, ii), nfrefs(ii)] = ref_select(Ytmp);
%         display(['!Done #', num2str(ii)])
    end
end





