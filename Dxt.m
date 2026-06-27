function d = Dxt(u)
%--------------------------------------------------------------------------
%   This program is computing the transpose of the gradient operator 
%   according to x-axis using backward difference.
%    
%   Usage: d = Dxt(u);
%
%   Inputs: 
%       - u: 2d data
%
%   Outputs: 
%       - d: transpose gradient u according to x-axis
% 
%   Code by: Xiaohao Cai
%   Last updated: 10/11/2012 
%--------------------------------------------------------------------------

[rows,cols] = size(u); 
d = zeros(rows,cols);
d(:,1:cols-1) = u(:,1:cols-1)-u(:,2:cols);

% use periodic boundary
d(:,cols) = u(:,cols)-u(:,1);

% %use reflective boundary
% d(:,cols)=0;
% d(:,cols)=d(:,cols-1);