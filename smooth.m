close all; clear all; clc;
addpath('data');

% load 4phase-data.mat
% Clean = new_Img;
% Img = imnoise(new_Img,'gaussian',0,0.25);
% Img=Img-min(Img(:));Img=Img/max(Img(:));
% k = 4;

load shape_multi_clean
Img=Img-min(Img(:));Img=Img/max(Img(:));
Clean = Img;
Img = double(rgb2gray(imread('given_new_2.png')));
% Img = imnoise(Img,'gaussian',0,0.1);
Img=Img-min(Img(:));Img=Img/max(Img(:));

k=4;
% 
iter_min = 2;
iter_max = 300;
[m,n] = size(Img);
p = m*n;H = 1;
Htg = H'*Img; HtH = H'*H;
t = 0.5;s = 1/16/t;
R = ClassFoImDiff('cir');
psf = [0 -1 0;-1 4 -1;0 -1 0];
Lap = ClassBlurMatrix(psf,[m,n]);
value_Lap = Lap*eye(m,n);
value_Lap(abs(value_Lap)<0.001) = 0;

%% initial values
% x0 = 0.5 * Img;
% xi = R * x0;
% du = R * x0;
Sigma = zeros(m,n);
alpha1 = 0.1;
alpha2 = 1;
% alpha1 = 0.01;
% alpha2 = 1;
sigma0 = 10;

lambda_array = zeros(iter_max,1);
mu_array = zeros(iter_max,1);
wx=0*Img;wy=wx;bx=wx;by=wx;
uker = 0*Img;
uker(1,1) = 4;uker(1,2)=-1;uker(2,1)=-1;
uker(end,1)=-1;uker(1,end)=-1;
x0 = 0.5 * Img;
%% update
for iter = 1:iter_max
    % update beta
    norma_beta = trace((Img-H*x0)*(Img-H*x0)');
    trazadiagonal = diag(sum(Sigma));
    trace_beta = trace(H'*H*trazadiagonal);
    beta = p/(norma_beta+trace_beta);
    
    % update lambda
    norma_lambda = sum(diag(value_Lap*(x0'*x0),0));
    trace_lambda = sum(diag(value_Lap * trazadiagonal,0));
    
    lambda = (2*alpha1*p)/(norma_lambda+trace_lambda);
    lambda_array(iter) = lambda;
    
    Du = R*x0;
    Dx = Du(:,:,:,1);Dy = Du(:,:,:,2);
    W_u = 1./(sqrt(Dx.^2+Dy.^2)+1e-3);
    W_u_m = mean(W_u(:));
%     wu_Lap = mean(W_u(:)) .* (value_Lap);
    
    
    % update mu
    trace_mu = sum(sum(Sigma .* (value_Lap) .* W_u));
    mu_traza = (Dx.^2+Dy.^2) .* W_u + abs(trace_mu)/p;
    tmp = mu_traza.^(0.5);
    mu = (p*alpha2)/(sum(tmp(:)));
    mu_array(iter) = mu;
    
    % update alpha1,alpha2
    alpha1 = (sum(lambda_array)/sum(lambda_array ~= 0)) * (norma_lambda+trace_lambda) / (2*p);
    alpha2 = (sum(mu_array)/sum(mu_array ~= 0))* (sum(tmp(:))) / p;
    
% %     % update the dual variable
%     xiold = xi;
%     z = xiold - s * du;
%     zz = max(sqrt(sum(sum(z.^2,4),3)),1);
%     xi = z./repmat(zz,[1 1 size(z,3) size(z,4)]);
%     
%     % update the primal variable
%     beta = 1;lambda = 3;mu = 5;
%     rhs = beta * t * Htg + x0 + mu * t * (R'*xi); 
%     A = beta * t *HtH + lambda * t * Lap + 1;
%     x = A\rhs;
%     
%     Sigma = 1./(beta * t *HtH+lambda * t * (value_Lap) + 1);
%     
%     % update the dual variable
%     du = R * x; 
%     z = xiold - s * du; 
%     zz = max(sqrt(sum(sum(z.^2,4),3)),1);
%     xi = z./repmat(zz,[1 1 size(z,3) size(z,4)]);

%     AS = beta * HtH + lambda * Lap + 1;
%     Sigma = AS \ eye(m,n);
%     Sigma = 1./(beta * t .* HtH+lambda * t .* value_Lap+1);
    
%       sigma = 2 / mu;
    
    sigma = mu * sigma0;   
    AS1 = beta * HtH + (lambda + sigma) * Lap;
    Sigma = AS1 \ (1 *eye(m,n));
        
    rhs = beta.*Htg+sigma*Dxt(wx-bx)+sigma*Dyt(wy-by);
    A = beta.*HtH+(lambda+sigma)*fft2(uker);
    x = ifft2(fft2(rhs)./ A);
    
    temp1=Dxf(x)+bx;temp2=Dyf(x)+by; 
    [wx,wy]=shrink2(temp1,temp2,1./sigma0);
    
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

[u,~] = CVMST_Seg(Img,1,2,2);
figure,imshow(Img);
figure,imshow(u);
figure,imshow(x);

th = ThdKmeans(x,k);
th1 = ThdKmeans(u,k);
% seg = SegResultShow(Img,x,th,k);
SA_value = new_SA(Clean,x,k,th);
SA_value1 = new_SA(Clean,u,k,th1);
disp(['our method is ' num2str(SA_value)]);
disp(['SaT method is ' num2str(SA_value1)]);