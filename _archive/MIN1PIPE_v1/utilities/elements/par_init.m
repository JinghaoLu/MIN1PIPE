function [P, options] = par_init(m)
% [P, options] = par_init initialization for CNMF
%   modified from E Pnevmatikakis
%   Jinghao Lu 05/20/2016

    [d1, d2, ~] = size(m, 'reg');                                % dimensions of dataset

    P = preprocess_data(m);
    
    options = CNMFSetParms(...                      
    'd1',d1,'d2',d2,...                         % dimensions of datasets
    'search_method','ellipse',...               % search locations when updating spatial components
    'deconv_method','constrained_foopsi',...    % activity deconvolution method
    'temporal_iter',2,...                       % number of block-coordinate descent steps 
    'fudge_factor',0.98);                       % bias correction for AR coefficients
end