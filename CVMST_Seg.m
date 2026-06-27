function [u,i]=CVMST_Seg(Img,lambda,mu,sigma)
%--------------------------------------------------------------------------
%   This program is using the Split-Bregman method to solve a convex variant
%   of the Mumford-Shah model proposed in our paper [1].
%    
%   Usage: [uu,te,er]=CVMST_Seg(Img,lambda,mu,sigma);
%
%   Inputs: 
%       - Img: the given image
%       - lambda: the fidelity parameter
%       - mu: the smooth term parameter
%       - sigma: the Split-Bregman parameter
%
%   Outputs: 
%       - u: the solution of our convex variant of the Mumford-Shah model
%
%
%   The original algorithm was proposed in the following paper, if you use 
%   this file or package for your work, please refer to the following paper:
% 
%   [1] Xiaohao Cai, Raymond H. Chan, and Tieyong Zeng,
%       A Two-Stage Image Segmentation Method Using a Convex Variant 
%       of the Mumford-Shah Model and Thresholding, accepted for 
%       publication by SIAM Journal on Imaging Sciences.
%
%   Code by: Xiaohao Cai
%   Last updated: 10/11/2012 
%--------------------------------------------------------------------------




%%
% initialization
u=Img;wx=0*Img;wy=wx;bx=wx;by=wx;

% Build Kernels: use the fft algorithm: (5 point stencil)
uker = 0*Img;
uker(1,1) = 4;uker(1,2)=-1;uker(2,1)=-1;
uker(end,1)=-1;uker(1,end)=-1;  

% compute the fft of the left hand side of (3.6) in our paper
uker = lambda+(mu+sigma)*fft2(uker);

%%
for i=1:300
    % compute the right hand side of (3.6) in our paper
    rhs = lambda.*Img+sigma*Dxt(wx-bx)+sigma*Dyt(wy-by);
        
    u0=u;
    
    % solve (3.6) in our paper
    u = ifft2(fft2(rhs)./uker);
    err=norm(u-u0,'fro')/norm(u,'fro');
    
%     if mod(i,10)==0
%         disp(['iterations: ' num2str(i) '!  ' 'error is:   ' num2str(err)]);
%     end
    
    % check the stopping criterion
    if err<10^(-5)
        break;
    end
    
    % solve (3.7) in our paper
    temp1=Dx(u)+bx;temp2=Dy(u)+by; 
    [wx,wy]=shrink2(temp1,temp2,1./sigma);
    
    % update (3.3) in our paper
    bx=bx+Dx(u)-wx;
    by=by+Dy(u)-wy;  
end

% disp(['All iterations: ' num2str(i)]);

% linear streach, i.e. (3.8) in our paper
min_u=min(min(u));
max_u=max(max(u));
u=(u-min_u)/(max_u-min_u);




