function requireR2015b

% Copyright 2015 The MathWorks, Inc.
% Steve Eddins

if verLessThan('matlab','R2015b')
   throwAsCaller(MException('se:ig:RequiresR2015b', ...
      'MATLAB R2015b is required.'));
end
