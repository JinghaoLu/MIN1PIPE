function [velocityField] = BCH(velocityField, updateVelocityField)
% BAKER_CAMPBELL_HAUSDORFF_FORMULA Accumulates two velocity fields using the Baker-Campbell-Hausdorff-approximation
% 
%   [velocityField] = baker_campbell_hausdorff_formula(...
%     velocityField, updateVelocityField, approximationMode)
%
%   INPUT ARGUMENTS
%   velocityField          - Accumulated velocity field
%   updateVelocityField    - Update velocity field
%   approximationMode      - Approximation to use ('A','B' or 'C')
%
%   OPTIONAL INPUT ARGUMENTS
%   N/A
%
%   OUTPUT ARGUMENTS
%   velocityField          - Accumulated velocity field

%   Copyright (c) 2012 Daniel Forsberg
%   danne.forsberg@outlook.com
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%   adapted by Jinghao Lu, 05/11/2017

    dims = ndims(velocityField{1});
    
    [v_xdx, v_xdy, v_ydx, v_ydy] = jacobian_matrix2d(velocityField);
    [u_xdx, u_xdy, u_ydx, u_ydy] = jacobian_matrix2d(updateVelocityField);
    jv = {v_xdx,v_xdy,v_ydx,v_ydy};
    ju = {u_xdx,u_xdy,u_ydx,u_ydy};
    firstExtraTerm = lie_bracket(velocityField, updateVelocityField, jv, ju);
    
    [u_xdx, u_xdy, u_ydx, u_ydy] = jacobian_matrix2d(firstExtraTerm);
    ju = {u_xdx,u_xdy,u_ydx,u_ydy};
    secondExtraTerm = lie_bracket(velocityField, firstExtraTerm, jv, ju);
    
    for k = 1: dims
        velocityField{k} = velocityField{k} + ...
            updateVelocityField{k} + 0.5 * firstExtraTerm{k} + ...
            1 / 12 * secondExtraTerm{k};
    end
end