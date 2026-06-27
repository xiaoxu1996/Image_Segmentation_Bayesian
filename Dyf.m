function d = Dyf(u)
%--------------------------------------------------------------------------
%   This program is computing the gradient operator according to y-axis 
%   using backward difference.
%    
%   Usage: d = Dy(u);
%
%   Inputs: 
%       - u: 2d data
%
%   Outputs: 
%       - d: gradient u according to y-axis
% 
%   Code by: Xiaohao Cai
%   Last updated: 10/11/2012 
%--------------------------------------------------------------------------

[rows,cols] = size(u); 
d = zeros(rows,cols);
d(2:rows,:) = u(2:rows,:)-u(1:rows-1,:);

% use periodic boundary
d(1,:) = u(1,:)-u(rows,:);

% %use reflective boundary
% d(1,:) = 0;
% d(1,:) = d(2,:);