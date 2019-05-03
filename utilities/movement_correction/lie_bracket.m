function [vectorField] = lie_bracket(vectorField,updateVectorField,jv,ju)
% LIE_BRACKET Computes the Lie bracket of vector fields
%
%   [vectorField] = lie_bracket(vectorField,updateVectorField)
%
%   INPUT ARGUMENTS
%   vectorField           - First vector field
%   updateVectorField     - Second vector field, assumed to be small
%
%   OPTIONAL INPUT ARGUMENTS
%   N/A
%
%   OUTPUT ARGUMENTS
%   vectorField 			- Final vector field

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

    v_xdx = jv{1};
    v_xdy = jv{2};
    v_ydx = jv{3};
    v_ydy = jv{4};
    u_xdx = ju{1};
    u_xdy = ju{2};
    u_ydx = ju{3};
    u_ydy = ju{4};
    temp{1} = (v_xdx .* updateVectorField{1} + v_xdy.*updateVectorField{2}) - ...
        (u_xdx .* vectorField{1} + u_xdy .* vectorField{2});
    temp{2} = (v_ydx .* updateVectorField{1} + v_ydy .* updateVectorField{2}) - ...
        (u_ydx .* vectorField{1} + u_ydy .* vectorField{2});
    vectorField = temp;
end