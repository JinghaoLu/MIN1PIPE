function [PNt, imout] = imregdft(imref, imcur)
% do phase correlation translation registration
%   Jinghao Lu, 10/20/2018

    [tf, tmp] = dftregistration(fft2(imref), fft2(imcur), 5);
    imout = ifft2(tmp);
    PNt = -tf(3: 4)';
end