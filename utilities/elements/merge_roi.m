function [roiout, sigout, seedsout, datasmthf, cutofff, pkcutofff] = merge_roi(m, roi, sig, seed, imax, datasmthf, cutofff, pkcutofff, ethres)
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
    roidenom = (repmat(sum(roig), nseed, 1) + repmat(sum(roig)', 1, nseed)) / 2;
    roicon = ((roig' * roig) ./ roidenom) > 0.5;
    roicon(1: nseed + 1: nseed ^ 2) = false;
    roicon = triu(roicon);
    
    %% compute cosine similarity (norm_inner) %%
    sigcorfn = sig_sim(datasmthf, cutofff, pkcutofff);
    sigcorfn = sigcorfn > ethres;
    
%     %%% salient peak keep unmerged %%%
%     [x, y] = ind2sub([d1, d2], seed);
%     salient_peak = false(nseed, nseed);
%     pthres = 0.01;
%     [X, Y] = meshgrid(1: d2, 1: d1);
%     for i = 1: nseed
%         xt = x(i);
%         yt = y(i);
%         for j = i + 1: nseed
%             yp = y(j);
%             xp = x(j);
%             d = sqrt((xt - xp) .^ 2 + (yt - yp) .^ 2);
%             if d < 2 * sz
%                 ltmp = interp2(X, Y, imax, linspace(yt, yp, 10), linspace(xt, xp, 10), 'linear', 5);
%                 ltmp = ltmp - linspace(ltmp(1), ltmp(end), 10);
%                 [mn, idt] = min(ltmp);
%                 if idt ~= 1 && idt ~= 10 && (ltmp(1) - mn) > pthres
%                     salient_peak(i, j) = true;
%                     salient_peak(j, i) = true;
%                 end
%             end
%         end
%     end

    %%% graph matrix finally used %%%
%     conmtx = roicon & sigcorfn & (~salient_peak);
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


