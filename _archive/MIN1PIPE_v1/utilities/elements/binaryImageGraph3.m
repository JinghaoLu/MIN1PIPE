function [g,nodenums] = binaryImageGraph3(bw,conn)
%binaryImageGraph3  Graph of foreground pixel connectivity in 3-D binary image
%
%   g = binaryImageGraph3(BW)
%   g = binaryImageGraph3(BW,conn)
%
%   [g,nodenums] = binaryImageGraph3(___)
%
%   g = binaryImageGraph3(BW) computes the 26-connected pixel
%   neighbor graph for the foreground pixels of a 3-D binary image.
%
%   g = binaryImageGraph3(BW,conn) computes the pixel neighbor graph
%   using the specified connectivity.
%
%   [g,nodenums] = binaryImageGraph3(___) also returns an array the
%   size of BW whose elements are the graph numbers corresponding to
%   each foreground pixel of BW.
%
%   EXAMPLES
%
%   Compute the pixel neighbor graph for the foreground pixels in a
%   3-D binary image.
%
%       bw = zeros(64,64,64);
%       bw(20:25,30:35,40:45) = 1;
%       g = binaryImageGraph3(bw)
%
%   Compute the 6-connected pixel neighbor graph for the foreground
%   pixels in a 3-D binary image.
%
%       bw = zeros(64,64,64);
%       bw(20:25,30:35,40:45) = 1;
%       g = binaryImageGraph3(bw,6)
%
%   MORE ABOUT CONNECTIVITY
%
%   3-D connectivity can be specified by the numbers 6, 18, or 26.
%   6-connected voxels share a common face. 18-connected voxels
%   share a common face or edge. 26-connected voxels share a common
%   face, edge, or vertex.
%
%   3-D connectivity can also be specified by a odd-size 3-D array
%   of 0s and 1s that is symmetric through its center element.
%
%   Note that a 3-D connectivity array can be larger than
%   3-by-3-by-3, which is a more general form than is allowed in
%   Image Processing Toolbox functions.
%
%   See also plotImageGraph, imageGraph, imageGraph3,
%   binaryImageGraph, adjacentRegionsGraph.

%   Copyright 2015 The MathWorks, Inc.
%   Steve Eddins

requireR2015b

% Input argument parsing and validation.
validateattributes(bw,{'numeric','logical'},{'3d'});
if ~islogical(bw)
   bw = (bw ~= 0);
end

default_conn = 26;
if nargin < 2
   conn = default_conn;
end
conn = checkConnectivity(conn);

% The connectivity array has already been verified to be symmetric
% about its center element. Because we are making an undirected
% graph, we don't need edges going in both directions, and so we
% discard (by setting to 0) half of the connectivity array.
conn(1:((end+1)/2)) = 0;

% Find the x, y, and z offsets (from center) of the nonzero elements
% in the connectivity array.
[conn_m,conn_n,conn_p] = size(conn);
[conn_offsets_y,conn_offsets_x,conn_offsets_z] = ind2sub([conn_m conn_n conn_p],find(conn));
conn_offsets_x = conn_offsets_x - ((conn_n + 1)/2);
conn_offsets_y = conn_offsets_y - ((conn_m + 1)/2);
conn_offsets_z = conn_offsets_z - ((conn_p + 1)/2);

% We have one node for each foreground pixels. Make a
% node_number_image for the purpose of mapping from an image pixel
% to the corresponding foreground pixel node number.
foreground_pixel_indices = find(bw);
nodenums = zeros(size(bw));
nodenums(bw) = (1:length(foreground_pixel_indices));

% Find the x, y, and z locations for each foreground pixel. We're
% going to store these locations as auxiliary information in the
% graph's node table.
[foreground_y,foreground_x,foreground_z] = ind2sub(size(bw),foreground_pixel_indices);

% Initialize an empty edge table.
edge_table = table(zeros(0,2),zeros(0,1),'VariableNames',{'EndNodes','Weight'});

% Build up the edge table by looking at one neighbor offset at a
% time.
for k = 1:length(conn_offsets_x)
   edge_table = [...
      edge_table
      binaryImageGraphNeighborEdges(bw,nodenums,...
      foreground_pixel_indices,foreground_x,foreground_y,foreground_z,conn_offsets_x(k),conn_offsets_y(k),conn_offsets_z(k))
      ];
end

node_table = table(foreground_x(:),foreground_y(:),foreground_z(:),...
   foreground_pixel_indices(:),'VariableNames',{'x','y','z','PixelIndex'});

g = graph(edge_table,node_table);

function edge_table = binaryImageGraphNeighborEdges(bw,node_number_image,foreground_pixel_indices,foreground_x,foreground_y,foreground_z,dx,dy,dz)

t = table;
t.PixelIndices = foreground_pixel_indices(:);
t.x = foreground_x(:);
t.y = foreground_y(:);
t.z = foreground_z(:);

% Calculate the neighbor locations.
neighbor_x = foreground_x + dx;
neighbor_y = foreground_y + dy;
neighbor_z = foreground_z + dz;

% Determine which neighbor candidates are out of bounds. Discard the
% corresponding entries in the table, t.
[bw_m,bw_n,bw_p] = size(bw);
out_of_bounds = (neighbor_x < 1) | (neighbor_x > bw_n) | ...
   (neighbor_y < 1) | (neighbor_y > bw_m) | ...
   (neighbor_z < 1) | (neighbor_z > bw_p);
t(out_of_bounds,:) = [];

% Determine which neighbor candidates belong to the image
% foreground.
neighbor_indices = sub2ind(size(bw),t.y+dy,t.x+dx,t.z+dz);
neighbor_not_foreground = ~bw(neighbor_indices);
t(neighbor_not_foreground,:) = [];
neighbor_indices(neighbor_not_foreground) = [];

% Return the edge table.
edge_table = table;
start = node_number_image(t.PixelIndices);
finish = node_number_image(neighbor_indices);
edge_table.EndNodes = [start(:) finish(:)];
edge_table.Weight = sqrt(dx^2 + dy^2 + dz^2) * ones(height(t),1);
