function [P,m] = preprocess_data(m,p,options)

% data pre-processing for:
% (i)   identifying and interpolating missing entries (assumed to have the
%       value NaN). Interpolated entries are passed back to Y.
% (ii)  identifying saturated pixels
% (iii) estimating noise level for every pixel
% (iv)  estimating global discrete time constants (if needed)
% This function replaces arpfit, present in the previous versions of the code.

% Author: Eftychios A. Pnevmatikakis
%           Simons Foundation, 2015

defoptions.noise_range = [0.25,0.5];            % frequency range over which to estimate the noise
defoptions.noise_method = 'logmexp';            % method for which to estimate the noise level
defoptions.block_size = [64,64];
defoptions.flag_g = false;                         % compute global AR coefficients
defoptions.lags = 5;                               % number of extra lags when computing the AR coefficients
defoptions.include_noise = 0;                      % include early lags when computing AR coefs
defoptions.split_data = 0;                         % split data into patches for memory reasons
defoptions.cluster_pixels = true;                  % cluster pixels into active or inactive

if nargin < 3 || isempty(options); options = defoptions; end
if nargin < 2 || isempty(p); p = 2; end
P.p = p;

if ~isfield(options,'noise_range'); options.noise_range = defoptions.noise_range; end
if ~isfield(options,'noise_method'); options.noise_method = defoptions.noise_method; end
if ~isfield(options,'block_size'); options.block_size = defoptions.block_size; end
if ~isfield(options,'flag_g'); options.flag_g = defoptions.flag_g; end
if ~isfield(options,'lags'); options.lags = defoptions.lags; end
if ~isfield(options,'include_noise'); options.include_noise = defoptions.include_noise; end; include_noise = options.include_noise;
if ~isfield(options,'split_data'); split_data = defoptions.split_data; else split_data = options.split_data; end
if ~isfield(options,'cluster_pixels'); cluster_pixels = defoptions.cluster_pixels; else cluster_pixels = options.cluster_pixels; end


%% interpolate missing data

[pixh, pixw, nf] = size(m, 'reg');
Y_interp = sparse(pixh * pixw, nf);
mis_data = [];
P.mis_values = full(Y_interp(mis_data));
P.mis_entries = mis_data;

%% estimate noise levels

fprintf('Estimating the noise power for each pixel from a simple PSD estimate...');
[sn, psx] = get_noise_fft(m, options);
P.sn = sn(:);
fprintf('  done \n');
end