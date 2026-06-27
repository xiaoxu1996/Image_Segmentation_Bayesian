close all; clear all; clc;
addpath('data');

load shape_multi_clean
Img=Img-min(Img(:));Img=Img/max(Img(:));
Clean = Img;
Img = double(rgb2gray(imread('given_new_2.png')));
% Img = imnoise(Img,'gaussian',0,0.1);
Img=Img-min(Img(:));Img=Img/max(Img(:));
k=4;

iter_min = 2;
iter_max = 300;
[m,n] = size(Img);
p = m*n;H = 1;
Htg = H'*Img; HtH = H'*H;

beta = 0.5;lambda = 2; 
sigma = 2;
mu = 1;
% sigma = mu * 10;
% sigma = 2;

num_samples = 1000;
Sigma = ones(m,n);


wx=0*Img;wy=wx;bx=wx;by=wx;
uker = 0*Img;
uker(1,1) = 4;uker(1,2)=-1;uker(2,1)=-1;
uker(end,1)=-1;uker(1,end)=-1;
x0 = 0.5 * Img;

R = ClassFoImDiff('cir');
psf = [0 -1 0;-1 4 -1;0 -1 0];
Lap = ClassBlurMatrix(psf,[m,n]);
value_Lap = Lap*eye(m,n);
value_Lap(abs(value_Lap)<0.001) = 0;

% 
% b1 = norm(H*Clean-Img,'fro');
% b2 = sum(diag(Clean * (R' * (R * Clean)),0));
% Du1 = R*Clean;
% Dx1 = Du1(:,:,:,1);Dy1 = Du1(:,:,:,2);
% b3 = norm(sqrt(Dx1.^2+Dy1.^2),'fro');
% prob1 = exp(-(beta/2) * b1 - (lambda/2) * b2 - mu * b3);

a1 = norm(H*x0-Img,'fro');
a2 = sum(diag(x0 * (R' * (R * x0)),0));
Du = R*x0;
Dx = Du(:,:,:,1);Dy = Du(:,:,:,2);
a3 = norm(sqrt(Dx.^2+Dy.^2),'fro');
prob = exp(-(beta/2) * a1 - (lambda/2) * a2 - mu * a3);

%% update
x_array = zeros(m,n);
xdata_array = zeros(m,n,num_samples);
nums = 0;
for t = 1:num_samples
%     xs = normrnd(x0,0.001);
    rhs = beta.*Htg+sigma*Dxt(wx-bx)+sigma*Dyt(wy-by);
    A = beta.*HtH+(lambda+sigma)*fft2(uker);
    xs = ifft2(fft2(rhs)./ A);
    b1 = norm(H*xs-Img,'fro');
    b2 = sum(diag(xs * (R' * (R * xs)),0));
    Du1 = R*xs;
    Dx1 = Du1(:,:,:,1);Dy1 = Du1(:,:,:,2);
    b3 = norm(sqrt(Dx1.^2+Dy1.^2),'fro');
    prob1 = exp(-(beta/2) * b1 - (lambda/2) * b2 - mu * b3);
    
    if prob1 > 100 * prob
        x_array = x_array + xs;
        nums = nums + 1;
        xdata_array(:,:,nums) = xs;
        prob = prob1;
    end
    
    temp1=Dxf(xs)+bx;temp2=Dyf(xs)+by; 
    [wx,wy]=shrink2(temp1,temp2,1./sigma);
    
    bx=bx+Dxf(xs)-wx;
    by=by+Dyf(xs)-wy;
end
mu_x = x_array / nums;
x = mu_x;
sigma_array = zeros(m,n);
for i = 1:nums
    xdata = xdata_array(:,:,nums);
    sigma_array = sigma_array + (xdata - x).^2;
end
sigma_array = sigma_array / nums;


% for iter = 1:iter_max
%     x_array = zeros(m,n);
%     for t = 1:num_samples
%         xs = normrnd(x0,(1./Sigma).^0.5);
%         a1 = norm(H*xs-Img,'fro');
%         a2 = sum(diag(xs * (R' * (R * xs)),0));
%         Du = R*xs;
%         Dx = Du(:,:,:,1);Dy = Du(:,:,:,2);
%         a3 = norm(sqrt(Dx.^2+Dy.^2),'fro');
%         prob = exp(-(beta/2) * a1 - (lambda/2) * a2 - mu * a3);
%         if prob > 1e-
% %         rhs = beta.*Htg+sigma*Dxt(wx-bx)+sigma*Dyt(wy-by)+Sigma * xs;
% %         A = beta.*HtH+(lambda+sigma)*fft2(uker) + fft2(Sigma);
% %         x = ifft2(fft2(rhs)./ A);
% %         Sigma = ifft2(fft2(eye(m,n))./ A);
% %         if (sum((x-xs).*(x-xs)) / sum(xs.*xs)) < 1e-5
% %             break;
% %         end
% %         xs = x;
%     end
% %     x = xs;
%     
% %     temp1=Dxf(x)+bx;temp2=Dyf(x)+by; 
% %     [wx,wy]=shrink2(temp1,temp2,mu./sigma);
% %     
% %     bx=bx+Dxf(x)-wx;
% %     by=by+Dyf(x)-wy;
% %     
% %     error = sum((x-x0).*(x-x0)) / sum(x0.*x0);
% %     if error < 1e-7 && iter>iter_min
% %         break; 
% %     end
% %     x0 = x;
%     
% end
%% Normalized
min_x=min(min(x));
max_x=max(max(x));
x=(x-min_x)/(max_x-min_x);


figure,imshow(x);

th = ThdKmeans(x,k);
% seg = SegResultShow(Img,x,th,k);
SA_value = new_SA(Clean,x,k,th);
disp(SA_value);





