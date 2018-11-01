function conn = checkConnectivity(conn)

% Copyright 2015 The MathWorks, Inc.
% Steve Eddins

if isequal(conn,4)
   conn = [0 1 0; 1 1 1; 0 1 0];
   
elseif isequal(conn,8)
   conn = ones(3,3);
   
elseif isequal(conn,6)
   conn = conndef(3,'minimal');
   
elseif isequal(conn,18)
   conn = cat(3,[0 1 0; 1 1 1; 0 1 0], ...
      ones(3,3), [0 1 0; 1 1 1; 0 1 0]);
   
elseif isequal(conn,26)
   conn = ones(3,3,3);
   
else
   % Assume that conn has been passed in array form and check its
   % validity.
   
   % Require that conn has odd dimensions, be 2-D or 3-D, and be
   % symmetric through the center element.
   
   if ndims(conn) > 3
      throwAsCaller(MException('se:ig:ConnectivityNumDims', ...
         'Connectivity must be 2-D or 3-D.'));
   end
   
   if any(mod(size(conn),2) == 0)
      throwAsCaller(MException('se:ig:ConnectivitySize',...
         'The connectivity size must be odd.'));
   end
   
   if ~isequal(conn(:), flipud(conn(:)))
      throwAsCaller(MException('se:ig:ConnectivitySymmetry',...
         'The connectivity must be symmetric through its center element.'));
   end
end