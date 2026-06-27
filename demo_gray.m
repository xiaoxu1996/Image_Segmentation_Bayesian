close all; clear all; clc;
addpath('data');

%% load 2Ps data
% Img=double(rgb2gray(imread('given_new_1.png')));
% Img=Img-min(Img(:));Img=Img/max(Img(:));
% 
% Clean = double(rgb2gray(imread('shape_clean.jpg')));
% Clean=Clean-min(Clean(:));Clean=Clean/max(Clean(:));
% k = 2;
% 
% th_Clean = ThdKmeans(Clean,k);
% Clean(Clean>=th_Clean) = 1;
% Clean(Clean< th_Clean) = 0;
% 
% 
% %% obtain the smoothed image and the parameters
% [x,beta,lambda,mu] = SaT_Bayesian_Seg(Img,Clean,k);
% 
% %% show the result and SA value
% th = ThdKmeans(x,k);
% seg = SegResultShow(Img,x,th,k);
% SA_value = new_SA(Clean,x,k,th);
% disp(SA_value);



%% load 4Ps data
load shape_multi_clean
Img=Img-min(Img(:));Img=Img/max(Img(:));
Clean = Img;
Img = double(rgb2gray(imread('given_new_2.png')));
Img=Img-min(Img(:));Img=Img/max(Img(:));
k=4;

%% obtain the smoothed image and the parameters
tic;
[x,beta,lambda,mu,para_time,alpha1_array,alpha2_array] = SaT_Bayesian_Seg(Img,Clean,4);
toc;


%% show the result and SA value
th = ThdKmeans(x,k);
seg = SegResultShow(Img,x,th,k);
SA_value = new_SA(Clean,x,k,th);
disp(SA_value);

%% data test parameter-nonnoise
% [Img,~] = imread('phase.bmp');
% Img = double(Img);
% 
% % Img=Img-min(Img(:));Img=Img/(max(Img(:))-min(Img(:)));
% load 3phase-lab.mat
% Clean = label;
% k = 3;

% [Img,~] = imread('CV8.bmp');
% Img = double(Img);
% load CV8phase-lab.mat
% Clean = label;
% k = 3;

% load 4phase-data-2.mat
% Clean = new_Img;
% Img = new_Img * 255;
% % Img = imnoise(new_Img,'gaussian',0,0.0001);
% % Img = (Img-min(Img(:)))/(max(Img(:))-min(Img(:)));
% k = 4; 

% load 5phase-lab.mat
% Clean = label;
% Img = imread('cam0836clean.bmp');
% 
% Img = rgb2gray(Img);
% Img = double(Img);
% k = 5;

% [x,alpha1_array,alpha2_array,SA_array] = SaT_Bayesian_Seg_array(Img,Clean,k);
% [x,beta,lambda,mu,para_time,alpha1_array,alpha2_array,SA_array] = SaT_Bayesian_Seg(Img,Clean,k);
% th = ThdKmeans(x,k);
% seg = SegResultShow(Img,x,th,k);
% SA_value = new_SA(Clean,x,k,th);
% disp(SA_value);