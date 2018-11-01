function I = iminterpolate(I, sx, sy)
% image interpolate for LogDemons
%   Jinghao Lu, 06/11/2016

    %%% Find update points on moving image %%%
    [x, y] = ndgrid(0: (size(I, 1) - 1), 0: (size(I, 2) - 1)); % coordinate image
    x_prime = x + sx; % updated x values (1st dim, rows)
    y_prime = y + sy; % updated y values (2nd dim, cols)
    
    %%% Interpolate updated image %%%
    I = interpn(x, y, I, x_prime, y_prime, 'linear', 0); % moving image intensities at updated points
end