function m = frame_stab(m)
% frame stabilization 
%   Jinghao Lu, 10/28/2018

    %%% initialize %%%
    [pixh, pixw, nf] = size(m, 'reg');
    sigmas = [1, 1, 1];
    hsizes = 2 * ceil(2 * sigmas) + 1;
    [hcol,hrow,hslc] = createSeparableGaussianKernel(sigmas, hsizes);
    ttype = class(m.reg(1, 1, 1));
    stype = parse_type(ttype);
    nsize = pixh * pixw * nf * stype * 10; %%% heuristic size of algorithm %%%
    
    %%% first spatial %%%
    nbatch = batch_compute(nsize);
    idbatch = round(linspace(0, nf, nbatch + 1));
    
    for i = 1: nbatch
        tmp = m.reg(1: pixh, 1: pixw, idbatch(i) + 1: idbatch(i + 1));
        tmp = convn(tmp, hcol, 'same');
        tmp = convn(tmp, hrow, 'same');
        m.reg(1: pixh, 1: pixw, idbatch(i) + 1: idbatch(i + 1)) = tmp;
    end
    
    %%% second temporal %%%
    nbatch = batch_compute(nsize);
    idbatch = round(linspace(0, pixh, nbatch + 1));
    
    for i = 1: nbatch
        tmp = m.reg(idbatch(i) + 1: idbatch(i + 1), 1: pixw, 1: nf);
        tmp = convn(tmp, hslc, 'same');
        m.reg(idbatch(i) + 1: idbatch(i + 1), 1: pixw, 1: nf) = tmp;
    end
end

function [hcol,hrow,hslc] = createSeparableGaussianKernel(sigma, hsize)

isIsotropic = all(sigma==sigma(1)) && all(hsize==hsize(1));

hcol = images.internal.createGaussianKernel(sigma(1), hsize(1));

if isIsotropic
    hrow = hcol;
    hslc = hcol;
else
    hrow = images.internal.createGaussianKernel(sigma(2), hsize(2));
    hslc = images.internal.createGaussianKernel(sigma(3), hsize(3));
end

hrow = reshape(hrow, 1, hsize(2));
hslc = reshape(hslc, 1, 1, hsize(3));

end
