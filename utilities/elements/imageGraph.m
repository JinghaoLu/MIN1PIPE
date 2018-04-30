function g = imageGraph(sz,conn)
%imageGraph  Graph of all image pixels
%
%   g = imageGraph(sz)
%   g = imageGraph(sz,conn)
%
%   g = imageGraph(sz) computes the 8-connected pixel neighbor graph for
%   a 2-D image with size specified by the two-element vector sz.
%
%   g = imageGraph(sz,conn) computes the pixel neighbor graph using the
%   specified connectivity.
%
%   EXAMPLES
%
%   Compute and plot the pixel neighbor graph for a 4-by-5 image.
%
%       g = imageGraph([4 5])
%       plotImageGraph(g)
%
%   Compute and plot the 4-connected pixel neighbor graph for a
%   4-by-5 image.
%
%       g4 = imageGraph([4 5],4)
%       plotImageGraph(g4)
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
%   See also plotImageGraph, imageGraph3, binaryImageGraph,
%   binaryImageGraph3, adjacentRegionsGraph.

%   Copyright 2015 The MathWorks, Inc.
%   Steve Eddins

requireR2015b

% Input argument parsing and validation.
validateattributes(sz,{'numeric'},{'row','integer','nonnegative',...
   'numel',2});

default_conn = 8;
if nargin < 2
   conn = default_conn;
end
conn = checkConnectivity(conn);

g = binaryImageGraph(true(sz),conn);
