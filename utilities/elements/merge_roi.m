function [roiout, sigout, seedsout, datasmthf, cutofff, pkcutofff] = merge_roi(m, roi, sig, seed, datasmthf, cutofff, pkcutofff, ethres)
% [roiout, sigout, seedsout, datasmthf, cutofff, pkcutofff] = merge_roi new merging criteria for
%   extracted rois
%   Jinghao Lu 06/10/2016

    hmerge = tic;
    %% initialization %%
    %%% initialize parameters %%%
    if nargin < 8 || isempty(ethres)
        defpar = default_parameters;
        ethres = defpar.merge_roi_corrthres;
    end
        
    nseed = length(seed);
    [d1, d2, d3] = size(m, 'reg');

    %% compute overlapping rois %%
    %%% this should be done in terms of gaussian smoothed rois %%%
    roig = roi_gauss(reshape(full(roi), d1, d2, nseed));
    roicon = (roig' * roig) > 0;
    roicon(1: nseed + 1: nseed ^ 2) = false;
    roicon = triu(roicon);
    
    %% compute cosine similarity (norm_inner) %%
    sigcorfn = sig_sim(datasmthf, cutofff, pkcutofff);
    sigcorfn = sigcorfn > ethres;
    
    %%% graph matrix finally used %%%
    conmtx = roicon & sigcorfn;
    
    %% get the graph %%
    iduse = merge_unit(conmtx, datasmthf);
    
    %% update variables %%
    roiout = roi(:, iduse);
    sigout = sig(iduse, :);
    seedsout = seed(iduse);
    datasmthf = datasmthf(iduse, :);
    cutofff = cutofff(iduse);
    pkcutofff = pkcutofff(iduse);
    time = toc(hmerge);
    disp(['Done merge roi, ', num2str(time), ' seconds'])
end


