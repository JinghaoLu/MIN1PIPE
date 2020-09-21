function [idusef, datasmthf, datusef, cutofff, pkcutofff] = seeds_merge(imax, iduse, datuse, datasmth, cutoff, pkcutoff, sz, corrthres)
% merge seeds within a same ROI
%   Jinghao Lu 05/20/2016

    %% initialization %%
    %%% initialize parameters %%%
    if nargin < 7 || isempty(sz)
        defpar = default_parameters;
        sz = defpar.neuron_size;
    end
    
    if nargin < 8 || isempty(corrthres)
        defpar = default_parameters;
        corrthres = defpar.pix_select_corrthres;
    end
    
    [pixh, pixw] = size(imax);
    nsd = length(iduse);
    
    %% compute signal similarity %%
    corrdata = sig_sim(datasmth, cutoff, pkcutoff, false);
    corrdatat = corrdata > corrthres; %%% empirical and liberal threshold %%%
    
    %% compute seeds distances/overlapping of potential neuron patches %%
    ovlpdata = zeros(size(corrdata));
    [x, y] = ind2sub([pixh, pixw], iduse);
    for i = 1: nsd
        xc = x(i);
        yc = y(i);
        ovlpdata(i, :) = sqrt((x - xc) .^ 2 + (y - yc) .^ 2)';
    end
    ovlpdata = triu(ovlpdata < 1.5 * sz) - diag(ones(1, nsd));    
        
    %%% salient peak keep unmerged %%%
    salient_peak = false(nsd, nsd);
    pthres = 0.005;
    [X, Y] = meshgrid(1: pixw, 1: pixh);
    for i = 1: nsd
        xt = x(i);
        yt = y(i);
        for j = i + 1: nsd
            yp = y(j);
            xp = x(j);
            d = sqrt((xt - xp) .^ 2 + (yt - yp) .^ 2);
            if d < 2 * sz
%                 ltmp = improfile(imax, [yt; yp], [xt; xp]);
                ltmp = interp2(X, Y, imax, linspace(xt, xp, 10), linspace(yt, yp, 10));
                ltmp = ltmp - linspace(ltmp(1), ltmp(end), 10);
                [mn, idt] = min(ltmp);
                if idt ~= 1 && idt ~= 10 && (ltmp(1) - mn) > pthres
                    salient_peak(i, j) = true;
                    salient_peak(j, i) = true;
                end
            end
        end
    end
    
    %%% graph matrix finally used %%%
    corrdatat = corrdatat & ovlpdata & (~salient_peak);
%     corrdatat = corrdatat & ovlpdata;
    
    %% get the graph %%
    iduset = seeds_merge_unit(corrdatat, imax, iduse);
    
    %% update variables %%
    idusef = iduse(iduset);
    datasmthf = datasmth(iduset, :);
    datusef = datuse(iduset, :);
    cutofff = cutoff(iduset);
    pkcutofff = pkcutoff(iduset);
end