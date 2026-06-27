function [x,beta,lambda,mu,para_time,alpha1_array,alpha2_array,SA_array] = SaT_Bayesian_Seg(Img,Clean,k)
% input: Img (Image to be segmented)
% output: x (the smoothed image)
%         beta (the beta parameter)
%         lambda (the lambda_1 parameter)
%         mu (the lambda_2 parameter)


%% variables set
iter_min = 2;
iter_max = 150;
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
x0 = 0.8 * Img;
xi = R * x0;
du = R * x0;
Sigma = zeros(m,n);
alpha1 = 0.1;
alpha2 = 0.1;

lambda_array = zeros(iter_max,1);
mu_array = zeros(iter_max,1);
mu_array(iter_max) = 0;
alpha1_array = zeros(iter_max,1);
alpha2_array = zeros(iter_max,1);
SA_array = zeros(iter_max,1); 
para_time = 0;
%% update
for iter = 1:iter_max
    
    alpha1_array(iter) = alpha1;
    alpha2_array(iter) = alpha2;
    
    % update beta    
    t00 = clock;
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
    
    % update mu
    trace_mu = sum(sum(Sigma .* (value_Lap) .* W_u));
    mu_traza = (Dx.^2+Dy.^2) .* W_u + abs(trace_mu)/p;
    tmp = mu_traza.^(0.5);
%     tmp = mu_traza;
    mu = (p*alpha2)/(sum(tmp(:)));
    mu_array(iter) = mu;
    
    % update alpha1,alpha2
    alpha1 = (sum(lambda_array)/sum(lambda_array ~= 0)) * (norma_lambda+trace_lambda) / (2*p);
    alpha2 = (sum(mu_array)/sum(mu_array ~= 0))* (sum(tmp(:))) / p;
    
    
    
    para_time = para_time + etime(clock,t00);
    
    % update the dual variable
    xiold = xi;
    z = xiold - s * du;
    zz = max(sqrt(sum(sum(z.^2,4),3)),1);
    xi = z./repmat(zz,[1 1 size(z,3) size(z,4)]);
    
    % update the primal variable
    rhs = beta * t * Htg + x0 + mu * t * (R'*xi); 
    A = beta * t *HtH + lambda * t * Lap + 1;
    x = A\rhs;
    
    Sigma = 1./(beta * t *HtH+lambda * t * (value_Lap) + 1);
    
    % update the dual variable
    du = R * x; 
    z = xiold - s * du; 
    zz = max(sqrt(sum(sum(z.^2,4),3)),1);
    xi = z./repmat(zz,[1 1 size(z,3) size(z,4)]);
    
    error = sum((x-x0).*(x-x0)) / sum(x0.*x0);
    if error < 1e-7 && iter>iter_min
        break; 
    end
    x0 = x;
    if Clean ~= 0
        min_x=min(min(x));
        max_x=max(max(x));
        Sa_x = (x-min_x)/(max_x-min_x);
        th = ThdKmeans(Sa_x,k);
        SA_array(iter) = new_SA(Clean,Sa_x,k,th);
    end  
end
disp(iter);
%% Normalized
min_x=min(min(x));
max_x=max(max(x));
x=(x-min_x)/(max_x-min_x);
