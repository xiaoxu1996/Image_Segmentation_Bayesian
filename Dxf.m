function d = Dxf(u)
%--------------------------------------------------------------------------
%   This program is computing the gradient operator according to x-axis 
%   using backward difference.
%    
%   Usage: d = Dx(u);
%
%   Inputs: 
%       - u: 2d data
%
%   Outputs: 
%       - d: gradient u according to x-axis
% 
%   Code by: Xiaohao Cai
%   Last updated: 10/11/2012 
%--------------------------------------------------------------------------


[rows,cols] = size(u);
d = zeros(rows,cols);
d(:,2:cols) = u(:,2:cols)-u(:,1:cols-1);
% use periodic boundary
d(:,1) = u(:,1)-u(:,cols);
% %use reflective boundary
% d(:,1) = 0;
% d(:,1) = d(:,2);