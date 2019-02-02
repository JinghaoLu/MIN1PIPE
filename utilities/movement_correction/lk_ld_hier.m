function [mxfn, xff, ldf, idf] = lk_ld_hier(mxall, mxalla, pixs, scl, sigma_x, sigma_f, sigma_d, sclld, maskc)
% LK-LogDemons hierarchical registration
%   Jinghao Lu, 02/27/2018

    if nargin < 2 || isempty(mxalla)
        mxalla = mxall;
    end
    
    if nargin < 3 || isempty(pixs)
        [pixh, pixw, ~] = size(mxall);
        pixs = min(pixh, pixw);
    end
    
    if nargin < 4 || isempty(scl)
        defpar = default_parameters;
        scl = defpar.mc_scl;
    end
    
    if nargin < 5 || isempty(sigma_x)
        defpar = default_parameters;
        sigma_x = defpar.mc_sigma_x;
    end
    
    if nargin < 6 || isempty(sigma_f)
        defpar = default_parameters;
        sigma_f = defpar.mc_sigma_f;
    end
    
    if nargin < 7 || isempty(sigma_d)
        defpar = default_parameters;
        sigma_d = defpar.mc_sigma_d;
    end
    
    if nargin < 8 || isempty(sclld)
        sclld = scl;
    end
    
    if nargin < 9 || isempty(maskc)
        maskc = true(size(mxall, 1), size(mxall, 2));
    end
    
    [pixh, pixw, ~] = size(mxall);
    
    %% LK first feedforward then full cluster %%

    %%% LK feedforward %%%
    [mxlkloop, xflkloop, idlkloop] = lk_loop(mxalla, pixs, scl, maskc);
    
    %%% LK full cluster %%%
    [mxlkclust, xflkclust, idlkclust] = lk_cluster(mxlkloop, pixs, scl, maskc);
    
    %%% combine the two LK layers %%%
    xflkcomb = [{xflkloop}, {xflkclust}];
    idlkcomb = [{idlkloop}, {idlkclust}];
    [xflkfnt, idlkfn] = lk_combine_layers(xflkcomb, idlkcomb);
    
    %%% convert the translation to affine form, and prepare sudo LogDemons transform %%%
    xflkfn = cell(1, length(xflkfnt));
    for i = 1: length(xflkfnt)
        xflkfn{i} = diag(ones(1, 3));
        nt = length(xflkfnt{i});
        for j = nt: -1: 1
            xflkfn{i} = xflkfnt{i}{j}.T * xflkfn{i};
        end
        xflkfn{i} = {affine2d(xflkfn{i})};
    end
    
    ldlkfn = cell(1, length(xflkfn));
    for i = 1: length(xflkfn)
        ldlkfn{i} = {{zeros(pixh, pixw, 2)}};
    end
    
    %%% prepare the combined images for the LogDemons layers %%%
    mx4ld = lk_warp_layers(squeeze(mat2cell(mxall, pixh, pixw, ones(1, size(mxall, 3)))), xflkfn);
    mx4ld = reshape(cell2mat(mx4ld(:)'), pixh, pixw, length(mx4ld));
    
    mxfn1 = zeros(pixh, pixw, length(idlkfn));
    for i = 1: length(idlkfn)
        mxfn1(:, :, i) = max(mx4ld(:, :, idlkfn{i}), [], 3);
    end
    mxfn1 = normalize(mxfn1);

    %% LogDemons first feedforward then full cluster %%
    
    %%% LogDemons feedforward %%%
    [mxout, xfuse, lduse, iduse] = logdemons_loop(mxfn1, pixs, scl, sigma_x, sigma_f, sigma_d, maskc);
    
    %%% LogDemons full cluster %%%
    [mxoutt, xfuset, lduset, iduset] = logdemons_cluster(mxout, pixs, sclld, sigma_x, sigma_f, sigma_d, maskc);
    
    %%% combine the two LogDemons layers %%%
    xfcomb = [{xfuse}, {xfuset}];
    ldcomb = [{lduse}, {lduset}];
    idcomb = [{iduse}, {iduset}];
    [xfldfn, ldfn, idldfn] = logdemons_combine_layers(xfcomb, ldcomb, idcomb);

    %% final combine LK and LogDemons transforms %%
    xfcomb = [{xflkfn}, {xfldfn}];
    ldcomb = [{ldlkfn}, {ldfn}];
    idcomb = [{idlkfn}, {idldfn}];
    [xff, ldf, idf] = logdemons_combine_layers(xfcomb, ldcomb, idcomb);
    mxfn = logdemons_warp_layers(squeeze(mat2cell(mxall, pixh, pixw, ones(1, size(mxall, 3)))), xff, ldf);
    mxfn = reshape(cell2mat(mxfn(:)'), pixh, pixw, length(mxfn));
end