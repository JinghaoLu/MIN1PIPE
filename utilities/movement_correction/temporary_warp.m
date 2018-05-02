function mxout = temporary_warp(xfallfn, idclustfn, mxintt)
% temporary warp original frames
%   Jinghao Lu, 10/11/2018

    %%% initialize %%%
    [pixh, pixw, ~] = size(mxintt);
    
    %%% combine layers %%%
    [xfuse, ~] = lk_combine_layers(xfallfn, idclustfn);

    %%% complete transform matrices %%%
    xflkfn = cell(1, length(xfuse));
    for i = 1: length(xfuse)
        xflkfn{i} = diag(ones(1, 3));
        nt = length(xfuse{i});
        for j = nt: -1: 1
            xflkfn{i} = xfuse{i}{j}.T * xflkfn{i};
        end
        xflkfn{i} = {affine2d(xflkfn{i})};
    end

    %%% final warp %%%
    mxout = lk_warp_layers(squeeze(mat2cell(mxintt, pixh, pixw, ones(1, size(mxintt, 3)))), xflkfn);
    mxout = reshape(cell2mat(mxout(:)'), pixh, pixw, length(mxout));
end
