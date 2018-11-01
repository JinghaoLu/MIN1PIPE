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
        
    %%% graph matrix finally used %%%
    corrdatat = corrdatat & ovlpdata;
    
    %% get the graph %%
    iduset = seeds_merge_unit(corrdatat, imax, iduse);
    
    %% update variables %%
    idusef = iduse(iduset);
    datasmthf = datasmth(iduset, :);
    datusef = datuse(iduset, :);
    cutofff = cutoff(iduset);
    pkcutofff = pkcutoff(iduset);
end