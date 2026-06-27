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

%% initial values
beta = 10;lambda = 90; 
mu = 2;
sigma = mu * 10;

wx=0*Img;wy=wx;bx=wx;by=wx;
uker = 0*Img;
uker(1,1) = 4;uker(1,2)=-1;uker(2,1)=-1;
uker(end,1)=-1;uker(1,end)=-1;
x0 = 0.5 * Img;
%% update
for iter = 1:iter_max
    
    rhs = beta.*Htg+sigma*Dxt(wx-bx)+sigma*Dyt(wy-by);
    A = beta.*HtH+(lambda+sigma)*fft2(uker);
    x = ifft2(fft2(rhs)./ A);
    
    temp1=Dxf(x)+bx;temp2=Dyf(x)+by; 
    [wx,wy]=shrink2(temp1,temp2,mu./sigma);
    
    bx=bx+Dxf(x)-wx;
    by=by+Dyf(x)-wy;
    
    error = sum((x-x0).*(x-x0)) / sum(x0.*x0);
    if error < 1e-7 && iter>iter_min
        break; 
    end
    x0 = x;
    
end
%% Normalized
min_x=min(min(x));
max_x=max(max(x));
x=(x-min_x)/(max_x-min_x);


figure,imshow(x);

th = ThdKmeans(x,k);
% seg = SegResultShow(Img,x,th,k);
SA_value = new_SA(Clean,x,k,th);
disp(SA_value);