function d = Dyt(u)
%--------------------------------------------------------------------------
%   This program is computing the transpose of the gradient operator 
%   according to y-axis using backward difference.
%    
%   Usage: d = Dyt(u);
%
%   Inputs: 
%       - u: 2d data
%
%   Outputs: 
%       - d: transpose gradient u according to y-axis
% 
%   Code by: Xiaohao Cai
%   Last updated: 10/11/2012 
%--------------------------------------------------------------------------

[rows,cols] = size(u); 
d = zeros(rows,cols);
d(1:rows-1,:) = u(1:rows-1,:)-u(2:rows,:);

% use periodic boundary
d(rows,:) = u(rows,:)-u(1,:);

% %use reflective boundary
% d(rows,:)=0;
% d(rows,:)=d(rows-1,:);