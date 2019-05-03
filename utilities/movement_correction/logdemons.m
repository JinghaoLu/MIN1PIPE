function [Mp, sx, sy] = logdemons(F, M, isgpu, nlevel, sigma_x, sigma_fluid, sigma_diffusion)
% Demons Registration
%
%   Simple matlab code for 2D image registration using the diffeomorphic log-demons algorithm 
%   Code is provided in order to help the understanding of the Demons algorithm
%
%   Original Algorithm in:
%   [1] - Symmetric Log-Domain Diffeomorphic Registration: A Demons-Based Approach
%         Vercauteren, Pennec, Perchant, Ayache -- MICCAI 2008, 754-761
%   [2] - Diffeomorphic demons: Efficient non-parametric image registration,
%         Vercauteren, Pennec, Perchant, Ayache -- NeuroImage 2009, (45)1:61-72
%
%   For a more recent work/survey, exploiting global shape characteristics 
%   (instead of the conventional local gradient-based approaches), consider citing 
%
%   [1] - Spectral Log-Demons: Diffeomorphic Image Registration with Very Large Deformations
%         Lombaert, Grady, Pennec, Ayache, Cheriet -- IJCV 2014, (107)3:254-271
%   Modified by Jinghao Lu, 06/16/2016

    %% Parameters
    if nargin < 4 || isempty(nlevel)
        nlevel = 3;   % multiresolution
    end
    
    if nargin < 5 || isempty(sigma_x)
        defpar = default_parameters;
        sigma_x = defpar.mc_sigma_x;
    end
    
    if nargin < 6 || isempty(sigma_fluid)
        defpar = default_parameters;
        sigma_fluid = defpar.mc_sigma_f;
    end
    
    if nargin < 7 || isempty(sigma_diffusion)
        defpar = default_parameters;
        sigma_diffusion = defpar.mc_sigma_d;
    end
    
    niter = 5;
    sigma_i = 1; % weight on similarity term
    do_display = 0;   % display iterations

    if isgpu && ~isa(F, 'gpuArray')
        FF = gpuArray(F);
        MM = gpuArray(M);
    else
        FF = F;
        MM = M;
    end

    if nlevel == 1
        %% Register
        opt = struct('niter', niter, 'sigma_fluid', sigma_fluid, 'sigma_diffusion', sigma_diffusion, 'sigma_i', sigma_i, 'sigma_x', sigma_x, 'do_display', do_display, 'do_plotenergy', 0);

        [Mp, sx, sy, ~, ~] = register(FF, MM, opt, isgpu);
    else
        %% Multiresolution
        vx = zeros(size(M)); % deformation field
        vy = zeros(size(M));
        for k=nlevel: -1: 1

            % downsample
            scale = 2 ^ -(k - 1);
            FFl = imresize(FF, scale);
            MMl = imresize(MM, scale);
            vxl = imresize(vx * scale, scale);
            vyl = imresize(vy * scale, scale);

            % register
            opt = struct('niter', niter,...
                'sigma_fluid', sigma_fluid,...
                'sigma_diffusion', sigma_diffusion,...
                'sigma_i', sigma_i,...
                'sigma_x', sigma_x,...
                'vx', vxl, 'vy', vyl,...
                'do_display', do_display, 'do_plotenergy', 1);
            [Mp, sxl, syl, vxl, vyl] = register(FFl, MMl, opt, isgpu);

            % upsample
            vx = imresize(vxl / scale, size(M));
            vy = imresize(vyl / scale, size(M));
        end
        sx = sxl;
        sy = syl;
    end
end

%% Register two images
function [Mp, sx, sy, vx, vy] = register(F, M, opt, isgpu)

    if nargin < 3
        opt = struct();  
    end
    if ~isfield(opt, 'sigma_fluid')  
        opt.sigma_fluid = 1.0;              
    end
    if ~isfield(opt, 'sigma_diffusion')
        opt.sigma_diffusion = 1.0; 
    end
    if ~isfield(opt, 'sigma_i')
        opt.sigma_i = 1.0;
    end
    if ~isfield(opt, 'sigma_x')
        opt.sigma_x = 1.0;
    end
    if ~isfield(opt, 'niter')
        opt.niter = 5;
    end
    if ~isfield(opt, 'stop_criterium')
        opt.stop_criterium = 0.01;
    end
    if ~isfield(opt, 'imagepad')
        opt.imagepad = 1.2;
    end
    if ~isfield(opt, 'vx')
        if isgpu
            opt.vx = zeros(size(M), 'gpuArray');
        else
            opt.vx = zeros(size(M));
        end
    end
    if ~isfield(opt, 'vy')
        if isgpu
            opt.vy = zeros(size(M), 'gpuArray');
        else
            opt.vy = zeros(size(M));
        end
    end
    
    %% padded image
    [F, ~] = imagepad(F, opt.imagepad, isgpu);
    [M, lim] = imagepad(M, opt.imagepad, isgpu);
    
    %% T is the deformation from M to F
    vx = imagepad(opt.vx, opt.imagepad, isgpu);
    vy = imagepad(opt.vy, opt.imagepad, isgpu);
    if isgpu
        e = zeros(1, opt.niter, 'gpuArray');
    else
        e = zeros(1, opt.niter);
    end
    e_min = 1e-0 + 100;      % Minimal energy
    
    %% Iterate update fields
    for iter = 1: opt.niter

        % Find update
        [ux, uy] = findupdate(F, M, vx, vy, opt.sigma_i, opt.sigma_x);

        % Regularize update
        ux = imgaussfilt(ux, opt.sigma_fluid);
        uy = imgaussfilt(uy, opt.sigma_fluid);
        
        % Compute step (e.g., max half a pixel)
        step = opt.sigma_x;

        %%% BCH %%%
        vfield = {vx, vy};
        ufield = {step * ux, step * uy};
        vcfield = BCH(vfield, ufield);
        vx = vcfield{1};
        vy = vcfield{2};
        
        % Regularize velocities
        vx = imgaussfilt(vx, opt.sigma_diffusion);
        vy = imgaussfilt(vy, opt.sigma_diffusion);
        
        % Get Transformation
        [sx, sy] = expfield(vx, vy);  % deformation field
        
        % Compute energy
        e(iter) = energy(F, M, sx, sy, opt.sigma_i, opt.sigma_x);
        if e(iter) < e_min
            sx_min = sx; sy_min = sy; % update best fields
            vx_min = vx; vy_min = vy; % update best fields
            e_min  = e(iter);
        end
        
        % Stop criterium
        if iter > 1 && abs(e(iter) - e(max(1, iter-5))) < e(1) * opt.stop_criterium
            break;
        end
    end
    
    %% Get Best Transformation
    if exist('vx_min', 'var')
        vx = vx_min;
        vy = vy_min;
        sx = sx_min;
        sy = sy_min;
    end
    
    %% Transform moving image
    Mp = iminterpolate(M, sx, sy);
    
    %% Unpad image
    Mp = Mp(lim(1): lim(2), lim(3): lim(4));
    vx = vx(lim(1): lim(2), lim(3): lim(4));
    vy = vy(lim(1): lim(2), lim(3): lim(4));
    sx = sx(lim(1): lim(2), lim(3): lim(4));
    sy = sy(lim(1): lim(2), lim(3): lim(4));
end

%% Find update between two images
function [ux,uy] = findupdate(F, M, vx, vy, sigma_i, sigma_x)

    % Get Transformation
    [sx, sy] = expfield(vx, vy);

    % Interpolate updated image
    M_prime = iminterpolate(M, sx, sy); % intensities at updated points
    
    % image difference
    diff = F - M_prime;
    
    % moving image gradient
    [gy,gx] = gradient_fast(M_prime);   % image gradient
    normg2 = gx .^ 2 + gy .^ 2;       % squared norm of gradient
%     area = size(M, 1) * size(M, 2); % area of moving image
    
    % update is Idiff / (||J||^2+(Idiff^2)/sigma_x^2) J, with Idiff = F(x)-M(x+s), and J = Grad(M(x+s));
    scale = diff ./ (normg2 + diff .^ 2 * sigma_i ^ 2 / sigma_x ^ 2);
    scale(normg2 == 0) = 0;
    scale(diff ==0) = 0;
    ux = gx .* scale;
    uy = gy .* scale;
    
    % Zero non overlapping areas
    ux(F == 0) = 0; 
    uy(F == 0) = 0;
    ux(M_prime == 0) = 0;
    uy(M_prime == 0) = 0;
end

%% Interpolate image
function I = iminterpolate(I, sx, sy)

    % Find update points on moving image
    [x, y] = ndgrid(0: (size(I, 1) - 1), 0: (size(I, 2) - 1)); % coordinate image
    x_prime = x + sx ; % updated x values (1st dim, rows)
    y_prime = y + sy ; % updated y values (2nd dim, cols)
    
    % Interpolate updated image
    I = interpn(x, y, I, x_prime, y_prime, 'linear', 0); % moving image intensities at updated points
end

%% Exponentiate vector field
function [vx, vy] = expfield(vx, vy)

    % Find n, scaling parameter
    normv2 = vx .^ 2 + vy .^ 2;
    m = sqrt(max(normv2(:)));
    n = ceil(log2(m / 0.5)); % n big enough so max(v * 2^-n) < 0.5 pixel)
    n = max(n, 0); % avoid null values
    
    % Scale it (so it's close to 0)
    vx = vx * 2 ^ -gather(n);
    vy = vy * 2 ^ -gather(n);

    % square it n times
    for i = 1: n
        [vx, vy] = compose(vx, vy, vx, vy);
    end
end

%% Compose two vector fields
function [vx, vy] = compose(ax, ay, bx, by)

    [x,y] = ndgrid(0: (size(ax, 1) - 1), 0: (size(ax, 2) - 1)); % coordinate image
    x_prime = x + ax; % updated x values
    y_prime = y + ay; % updated y values
    
    % Interpolate vector field b at position brought by vector field a
    bxp = interpn(x, y, bx, x_prime, y_prime, 'linear', 0); % interpolated bx values at x+a(x)
    byp = interpn(x, y, by, x_prime, y_prime, 'linear', 0); % interpolated bx values at x+a(x)

    % Compose
    vx = ax + bxp;
    vy = ay + byp;
end

%% Jacobian
function det_J = jacobian(sx, sy)

    % Gradients
    [gx_y, gx_x] = gradient_fast(sx);
    [gy_y, gy_x] = gradient_fast(sy);
    
    % Add identity
    gx_x = gx_x + 1; % zero displacement should yield a transformation T = Identity (points keep their positions)
    gy_y = gy_y + 1; % adding identity matrix here
    
    % Determinant
    det_J = gx_x .* gy_y - gy_x .* gx_y;
end


%% Get energy
function e = energy(F, M, sx, sy, sigma_i, sigma_x)

    % Intensity difference
    Mp = iminterpolate(M, sx, sy);
    diff2 = (F - Mp) .^ 2;
    area = size(M, 1) * size(M, 2);
    
    % Transformation Gradient
    jac = jacobian(sx, sy);
    
    % Three energy components
    e_sim = sum(diff2(:)) / area;
    %e_dist = sum((cx(:)-sx(:)).^2 + (cy(:)-sy(:)).^2) / area;
    e_reg = sum(jac(:) .^ 2) / area;
    
    % Total energy
    e = e_sim + (sigma_i ^ 2 / sigma_x ^ 2) * e_reg;
end

%% Pad image
function [I,lim] = imagepad(I, scale, isgpu)

    if nargin < 2
        scale = 2;
    end % default, pad image twice as big
    
    if isgpu
        Ip = zeros(ceil(size(I) * scale), 'gpuArray');
    else
        Ip = zeros(ceil(size(I) * scale));
    end
    lim = bsxfun(@plus, floor(size(I) * (scale - 1) / 2), [[1 1]; size(I)]); % image limits
    Ip(lim(1): lim(2), lim(3): lim(4)) = I;                              % padded image
    I = Ip;
end


