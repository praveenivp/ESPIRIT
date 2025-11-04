%% ESPIRiT Maps Demo
% This is a demo on how to generate 3D ESPIRiT maps. It is based on the paper
% Uecker et. al, MRM 2013 DOI 10.1002/mrm.24751 and BART. We control

%%
% Set parameters
addpath(genpath('C:\Users\pvalsala\Documents\Packages2\ESPIRiT'))
addpath(genpath('C:\Users\pvalsala\Documents\Packages2\arrShow'))
load brain_8ch.mat
calib_data=reshape(crop(single(DATA),[24,24,8]),[24 24 1 8]);
%calib_data is of size [Sx,Sy,Sz,Nc]
imSize=[200 200 1];
ksize = [1,1,1]*6; % kernel size
% Threshold for picking singular vercors of the calibration matrix
% (relative to largest singlular value.
eigThresh_1 = 0.01;
% threshold of eigen vector decomposition in image space. 
eigThresh_2 = 0.95*0; %set to 0 if you want to see null space!

%% Compute ESPIRiT EigenVectors
% Here we perform calibration in k-space followed by an eigen-decomposition
% in image space to produce the EigenMaps. 
% compute Calibration matrix, perform 1st SVD and convert singular vectors
% into k-space kernels
kernels=dat2Kernel3D(calib_data,ksize,eigThresh_1);

%%
nMaps=8;
tic
[EigenVecs, EigenVals] = kernelEig3D(kernels, imSize,nMaps,eigThresh_2);
toc
EigenVecs=flip(flip(flip(EigenVecs,1),2),5);

%%
% show eigen-values and eigen-vectors. The last set of eigen-vectors
% corresponding to eigen-values 1 look like sensitivity maps

figure, imshow3(abs(EigenVals),[],[1,nMaps]); 
title('Eigen Values in Image space');
colormap((gray(256))); colorbar;

figure, imshow3(abs(EigenVecs),[],[8,nMaps]); 
title('Magnitude of Eigen Vectors');
colormap(gray(256)); colorbar;

figure, imshow3(angle(EigenVecs),[],[8,nMaps]); 
title('Phase of Eigen Vectors');
colormap(jet(256)); colorbar;

%%
% project onto the eigenvectors. This shows that all the signal energy
% lives in the subspace spanned by the eigenvectors with eigenvalue 1.
% These look like sensitivity maps. 

im_coils =permute( ifft2c(DATA), [1 2 4 3]);
[ip,proj,null]=espirit_proj(im_coils,EigenVecs);


figure, subplot(311),imshow3(abs(ip),[],[1,8]); 
title('inner product');
colormap(sqrt(gray(256))); colorbar;

subplot(312),imshow3(abs(proj),[],[1,8]); 
title('projection');
colormap(sqrt(gray(256))); colorbar;

subplot(313),imshow3(abs(null),[],[1,8]); 
title(sprintf('null space for eig-threshold1=%.2f and %d kernel',eigThresh_1,size(kernels,5)));
colormap(sqrt(gray(256))); colorbar;

%% let's try full 3D
calib_data=readcfl('calib24');

imSize=[200 200 40];
ksize = [1,1,1]*6; % kernel size
eigThresh_1 = 0.01; 
eigThresh_2 = 0.95*0; %set to 0 if you want to see null space!
nMaps=8;

tic
kernels=dat2Kernel3D(calib_data,ksize,eigThresh_1);
toc

tic
[EigenVecs, EigenVals] = kernelEig3D(kernels, imSize,nMaps,eigThresh_2);
toc
EigenVecs=flip(flip(flip(flip(EigenVecs,1),2),3),50);


%%
function [ip,proj,null]=espirit_proj(im_coils,maps)

% innner product: what each maps 
ip=sum(im_coils.*conj(maps),4);

%proj:
proj=sum(ip.*maps,5);

null=im_coils-proj;
end
