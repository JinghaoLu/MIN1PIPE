function [m_out, pixh, pixw] = downsamp_unit(m_in, spatialr, ttype)
    if nargin < 6 || isempty(ttype)
        defpar = default_parameters;
        ttype = defpar.ttype;
    end
    [pixh, pixw, nf] = size(m_in, 'frame_allt');
    
    filename = m_in.Properties.Source;
    filename = [filename(1: end - 5), filename(end - 3: end)];
    m_out = matfile(filename, 'Writable', true);
    
    %% batch configuration %%
    stype = parse_type(ttype);
    nsize = pixh * pixw * nf * stype; %%% size of single %%%
    nbatch = batch_compute(nsize);
    ebatch = ceil(nf / nbatch);

    %%% extract batch-wise frame info %%%
    idbatch = [1: ebatch: nf, nf + 1];

    %% downsampling %%
    if isempty(gcp('nocreate'))
        parpool(feature('numCores'));
    end
    
    pixht = round(pixh * spatialr);
    pixwt = round(pixw * spatialr);
    for ib = 1: nbatch
        frame_all = m_in.frame_allt(1: pixh, 1: pixw, idbatch(ib): idbatch(ib + 1) - 1);
        frame_allt = zeros(pixht, pixwt, size(frame_all, 3), class(frame_all));
        parfor j = 1: size(frame_all, 3)
            frame_allt(:, :, j) = imresize(frame_all(:, :, j), [pixht, pixwt]);
        end
        frame_all = frame_allt;
        savef(filename, 2, 'frame_all')
    end
    pixh = pixht;
    pixw = pixwt;
end