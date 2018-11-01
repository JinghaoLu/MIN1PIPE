function [g,nodenums] = binaryImageGraph(bw,conn)
%binaryImageGraph  Graph of foreground pixel connectivity in binary image
%
%   g = binaryImageGraph(BW)
%   g = binaryImageGraph(BW,conn)
%
%   [g,nodenums] = binaryImageGraph(___)
%
%   g = binaryImageGraph(BW) computes the 8-connected pixel neighbor
%   graph for the foreground pixels of a 2-D binary image.
%
%   g = binaryImageGraph(BW,conn) computes the pixel neighbor graph
%   using the specified connectivity.
%
%   [g,nodenums] = binaryImageGraph(___) also returns a matrix the
%   size of BW whose elements are the graph numbers corresponding to
%   each foreground pixel of BW.
%
%   EXAMPLES
%
%   Compute and plot the pixel neighbor graph for the foreground
%   pixels in the image text.png (a sample image shipped with the
%   Image Processing Toolbox).
%
%       bw = imread('text.png');
%       g = binaryImageGraph(bw);
%       plotImageGraph(g)
%       % Zoom in
%       axis([60 85 30 45])
%
%   Compute and plot the 4-connected pixel neighbor graph for the
%   foreground pixels in the image text.png (a sample image shipped
%   with the Image Processing Toolbox).
%
%       bw = imread('text.png');
%       g = binaryImageGraph(bw,4);
%       plotImageGraph(g)
%       % Zoom in
%       axis([60 85 30 45])
%
%   MORE ABOUT CONNECTIVITY
%
%   Two-dimensional connectivity can be specified by the number 4 or
%   the number 8, specifying the traditional definitions of
%   4-connected or 8-connected pixel neighbors.
%
%   Two-dimensional connectivity can also be specified by a odd-size
%   matrix of 0s and 1s that is symmetric through its center
%   element. For example, the following connectivity matrix says
%   that a pixel is a connected neighbor of the pixel just above it
%   and the pixel just below it:
%
%       conn = [ ...
%                0  1  0
%                0  1  0
%                0  1  0  ]
%
%   And the following connectivity matrix says that a pixel is
%   a connected neighbor of the pixels to the upper right and lower
%   left.
%
%       conn = [ ...
%                0  0  1
%                0  1  0
%                1  0  0  ]
%
%   Note that the connectivity matrix can be larger than 3-by-3,
%   which is a more general form than is allowed in Image Processing
%   Toolbox functions.
%
%   See also plotImageGraph, imageGraph, imageGraph3,
%   binaryImageGraph3, adjacentRegionsGraph.

%   Copyright 2015 The MathWorks, Inc.
%   Steve Eddins

requireR2015b

% Input argument parsing and validation.
validateattributes(bw,{'numeric','logical'},{'2d'});
if ~islogical(bw)
   bw = (bw ~= 0);
end

default_conn = 8;
if nargin < 2
   conn = default_conn;
end
conn = checkConnectivity(conn);

[g,nodenums] = binaryImageGraph3(bw,conn);
g.Nodes.z = [];
