%% see coil profiles of compressed coil
reco_obj=recoVBVD('S:\KYB\AGKS\pvalsala\20210114_SpiralAcceltest_b0\meas_MID43_B0mapping_2mm_6E_BW250_FID22853.dat');


%%
ref=reco_obj.twix.refscan{''};

imsize=[reco_obj.NCol reco_obj.NLin reco_obj.NPar ];

t=(fftshift(fftshift(fft(fft(fft(ref,imsize(1),1),imsize(2),3),imsize(3),4),1),3));
as(t)

%%
% [fn,path]=uigetfile('S:\KYB\AGKS\pvalsala\20210114_SpiralAcceltest_b0\*.dat');
sp=SpiralReco('S:\KYB\AGKS\pvalsala\20210114_SpiralAcceltest_b0\meas_MID56_spiralout_i64_FID22866.dat');
figure,montage(permute(sp.coilSens,[2 3 4 1]));


[imc1,coilsens1]=adaptiveCombine2(reco_obj.img(:,:,:,:,1,1,1,1));

%%
im1=makemontage(permute(sp.coilSens,[2 3 1]));
%%
im2=makemontage(permute(coilsens1(:,:,:,9),[2 3 4 1]));




%%
%% ESPIRiT Maps Demo
% This is a demo on how to generate ESPIRiT maps. It is based on the paper
% Uecker et. al, MRM 2013 DOI 10.1002/mrm.24751. ESPIRiT is a method that
% finds the subspace of multi-coil data from a calibration region in
% k-space using a series of eigen-value decompositions in k-space and image
% space. 

ref_parfft=fftshift(fft(ref,[],4),4);
DATA=ref_parfft(:,:,:,8,1);
DATA=permute(DATA,[1 3 2]);

%%
[sx,sy,Nc] = size(DATA);
ncalib = 24; % use 24 calibration lines to compute compression
ksize = [6,6]; % kernel size


% Threshold for picking singular vercors of the calibration matrix
% (relative to largest singlular value.

eigThresh_1 = 0.02;

% threshold of eigen vector decomposition in image space. 
eigThresh_2 = 0.95;

% crop a calibration area
calib = crop(DATA,[ncalib,ncalib,Nc]);

%%
% Display coil images: 
im = ifft2c(DATA);

figure, imshow3(abs(im),[],[1,Nc]); 
title('magnitude of physical coil images');
colormap((gray(256))); colorbar;

figure, imshow3(angle(im),[],[1,Nc]); 
title('phase of physical coil images');
colormap('default'); colorbar;

%% Compute ESPIRiT EigenVectors
% Here we perform calibration in k-space followed by an eigen-decomposition
% in image space to produce the EigenMaps. 


% compute Calibration matrix, perform 1st SVD and convert singular vectors
% into k-space kernels

[k,S] = dat2Kernel(calib,ksize);
idx = max(find(S >= S(1)*eigThresh_1));

%% 
% This shows that the calibration matrix has a null space as shown in the
% paper. 

kdisp = reshape(k,[ksize(1)*ksize(2)*Nc size(k,4)]);
figure, subplot(211), plot([1:size(k,4)],S,'LineWidth',2);
hold on, 
plot([1:size(k,4)],S(1)*eigThresh_1,'r-','LineWidth',2);
plot([idx,idx],[0,S(1)],'g--','LineWidth',2)
legend('signular vector value','threshold')
title('Singular Vectors')
subplot(212), imagesc(abs(kdisp)), colormap(gray(256));
xlabel('Singular value #');
title('Singular vectors')


%%
% crop kernels and compute eigen-value decomposition in image space to get
% maps
[M,W] = kernelEig(k(:,:,:,1:idx),[100 100]);

%%
% show eigen-values and eigen-vectors. The last set of eigen-vectors
% corresponding to eigen-values 1 look like sensitivity maps

figure, imshow3(abs(W),[],[1,Nc]); 
title('Eigen Values in Image space');
colormap((gray(256))); colorbar;

figure, imshow3(abs(M),[],[Nc,Nc]); 
title('Magnitude of Eigen Vectors');
colormap(gray(256)); colorbar;

figure, imshow3(angle(M),[],[8,8]); 
title('Magnitude of Eigen Vectors');
colormap(jet(256)); colorbar;


%%
% project onto the eigenvectors. This shows that all the signal energy
% lives in the subspace spanned by the eigenvectors with eigenvalue 1.
% These look like sensitivity maps. 


% alternate way to compute projection is:
% ESP = ESPIRiT(M);
% P = ESP'*im;

P = sum(repmat(im,[1,1,1,Nc]).*conj(M),3);
figure, imshow3(abs(P),[],[1,Nc]); 
title('Magnitude of Eigen Vectors');
colormap(sqrt(gray(256))); colorbar;

figure, imshow3(angle(P),[],[1,Nc]); 
title('Magnitude of Eigen Vectors');
colormap((jet(256))); colorbar;



%%
% crop sensitivity maps 
maps = M(:,:,:,end).*repmat(W(:,:,end)>eigThresh_2,[1,1,Nc]);

figure, imshow3(abs(maps),[],[1,8]); 
title('Absolute sensitivity maps');
colormap((gray(256))); colorbar;

figure, imshow3(angle (maps),[],[1,8]); 
title('Phase of sensitivity maps');
colormap((jet(256))); colorbar;
%%
function im_montage= makemontage(im,sz)
im=reshape(im,size(im,1),size(im,2),[]);

if( ~exist('sz','var')||isempty(sz))
    sz=ones([2 1]);
    sz(1)=ceil(sqrt(size(im,3)));
    sz(2)=ceil(size(im,3)/sz(1));
end

im=padarray(im,[0 0 prod(sz)-size(im,3)],0,'post');
im_montage=reshape(permute(reshape(im,[size(im,1) size(im,2) sz(1) sz(2)]),[1 3 2 4]),[size(im,1)*sz(1) size(im,2)*sz(2)]);
if(isreal(im_montage))
figure,imagesc((im_montage)),colorbar,colomap('parula')
else
   figure,
   ax=subplot(121);
   imagesc(abs(im_montage)),colorbar,colormap(ax,'parula')
   ax=subplot(122);
   imagesc(angle(im_montage)),colorbar,colormap(ax,'hsv')
end

end

        function [cdata,V]=performCoilCompression(data)
          
                
                data=data(:,:).';
                [~,sval,vec]=svd(data,'econ');
        
                    tol = 0.95;
               
                NcCha = find((cumsum(diag(sval))/sum(diag(sval))) >tol,1,'first');
                V=vec(1:NcCha,:);
                fprintf('Compressing %d coils into %d coils with %d%% signal energy\n',size(sval,1),NcCha,tol*100);
                sz    = size(data);
                cdata   = V*data(:,:);
                cdata   = reshape(cdata,[NcCha sz(2:end)]);
                
        end