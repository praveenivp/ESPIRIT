%% Load file and caculate trajectories
% filename='S:\KYB\AGKS\pvalsala\Phantom\build_pencspiral_13112020\meas_MID98_peNc_spiralout_i64_msint_rfsp_R1_FID16643.dat';
% filename='S:\KYB\AGKS\pvalsala\Phantom\build_pencspiral_13112020\meas_MID102_peNc_spiralout_i64_rfsp_PAT3_A32_FID16647.dat';
% filename='S:\KYB\AGKS\pvalsala\Phantom\build_pencspiral_13112020\meas_MID100_peNc_spiralout_i64_rfsp_PAT2_A16_FID16645.dat';
[fn,path]=uigetfile('S:\KYB\AGKS\pvalsala\20210114_SpiralAcceltest_b0\*.dat');
filename=fullfile(path,fn);
twix=mapVBVD(filename);
parameter=getSpiralPara(twix);
[G_PRS,G_xyz,grad_MOM]=GetGradients(twix);
load('D:\Data\GIRF_fieldcamera\20200828_GIRFnew\PSF_time.mat')

%read phase slice to Phase read slice
G_corr=(GIRF_Correction(G_xyz,PSF_time,'isCrossTermPSFCorr',true));
G_PRS_corr=GradientXYZ2PRS(G_corr(:,2:4,:),twix);
parameter.grad_MOM=grad_MOM;
k=Grad2Traj(G_PRS_corr,parameter,'PE');
w=jacksonDCF2(k,parameter);
kmax=2*pi*(0.5/(parameter.Resolution*1e-3));
k_scaled=k./(2*kmax);
N=parameter.FOV(1)/parameter.Resolution;
N = [N,N];                  % size of the target image
GFFT3 = NUFFT(k_scaled(:),w(:),1, [0,0] , N,2);

%% Looping starts
para1=[3 5 7 9 11 13 15];
para2=[10 20 30 40 50 60 70 80];
% para1=[3 5 ];
% para2=[10 ];
im_spirit=zeros([N length(para1) length(para2)]);
for ii=1:length(para1)
for jj=1:length(para2)
%% SpiriT parameters
% kSize = [7,7];                  % SPIRiT Kernel size
% CalibSize = [80,80];            % size of the calibration region
kSize = [1,1]*para1(ii);                  % SPIRiT Kernel size
CalibSize = [1,1]*para2(jj);            % size of the calibration region

nIterCG = 30;                   % number of reconstruction iterations
CalibTyk = 0.02;                % Tykhonov regularization for calibration 0.01-0.05 recommended
lambda = 1;                     % Ratio between data and calibration consistency. 1 recommended when density compensation is used.
accel = 2;%parameter.PATmode;                      % Acceleration factor

settings= sprintf('R=%d  |   Kernel size =[%d %d] |  ACS size=[%d %d]\n' ,accel,kSize,CalibSize);
settings=strcat(settings,sprintf('nIterCG =%d | lambda =%d | ',nIterCG,lambda));

eigThresh_k = 0.02; % threshold of eigenvectors in k-space
eigThresh_im = 0.95; % threshold of eigenvectors in image space

%% load data and perform FOV shift
acs=twix.refscan{''};
pos_PRS=GradientXYZ2PRS(1e-3*[1 -1 -1].*parameter.slice{1}.Position,twix); %only work for head first-supine
B0_mod=exp(-1i*(real(k).*pos_PRS(1)+imag(k).*pos_PRS(2)));
% negative displacement added to real part of kTRaj moves up down in array show

% coil compression
acs=acs(:,:,:,end);
% [acs_comp,nCoil,V]= CoilcompressSVD(acs,[],0.01);
 acs_comp=permute(acs,[1 3 2]); % COL X LIN X X COIL
 nCoil=size(acs_comp,3);
acs_comp=bsxfun(@times, acs_comp,sqrt(w).*B0_mod);
%%
im = GFFT3'*double(acs_comp);
kData = fft2c(im);


% kernel = zeros([kSize,nCoil,nCoil]);
kCalib = crop(kData,[CalibSize,nCoil]); % crop center k-space for calibration

% prepare calibration matrix for calibration ( Phil Beatty calls this
%correlation values matrix. See thesis by Phil Beatty.)
%[AtA,] = corrMatrix(kCalib,kSize);
%for n=1:nCoil
%    disp(sprintf('Calibrating coil %d',n));
%	kernel(:,:,:,n) = calibrate(AtA,kSize,nCoil,n,CalibTyk);
%end
tic
kernel = calibSPIRiT(kCalib, kSize, nCoil, CalibTyk);
GOP = SPIRiT(kernel, 'image',N); % init the SPIRiT Operator

fprintf('Done Calibrating in %3.2f seconds\n',toc);

kernel_time=toc;
%% load undersampoled data

GFFT3 = NUFFT(k_scaled,w,1, [0,0] , N,2);
data=twix.image{''};
idx =(1:accel:size(k_scaled,2));
k_u = k_scaled(:,idx);
w_u = w(:,idx);  % use w_u = w(:,idx)*0+1; if you don't want to use density weighting
disp('generating nufft object for recon')
GFFT_u = NUFFT(k_u,w_u,1, [0,0], N,2);

%No compression
kData_u=permute(data(:,:,idx),[1 3 2]); % COL X LIN X X COIL
nCoil=size(kData_u,3);
%compress coils
% [kData_u,nCoil]= CoilcompressSVD(data(:,:,idx),V);
kData_u=double(bsxfun(@times, kData_u,sqrt(w_u).*B0_mod(:,idx)));

%% Undersample the data and prepare new NUFFT operator for it


% this may improve noise, but will converge slower. use
% larger lambda.
%  kData_u = double(ref1(:,idx,:)); % if using coilcompress and kernel has
%  same coil

im_zp = (GFFT_u'*kData_u)*accel;
res = im_zp;

disp('initiating reconstruction')
tic
res = cgNUSPIRiT(kData_u,res,GFFT_u,GOP,nIterCG*2,lambda);
fprintf('done SpiriT reconstructing in %3.2f seconds!\n',toc);
spirit_time=toc;


settings=strcat(settings,sprintf('kernel Time =%1.2f seconds | Spirit time =%1.2f seconds | ',kernel_time,spirit_time));
all_set{ii,jj}=settings;
% 
im_acs=sos(im);
im_zp=sos(im_zp);
im_spirit(:,:,ii,jj)=sos(res);
h=figure;
imshow(cat(2,im_acs/max(im_acs(:)),im_zp./max((im_zp(:))),im_spirit(:,:,ii,jj)./max(col(im_spirit(:,:,ii,jj)))),[]),

annotation(h,'textbox',...
    [0.365583333333333 0.217913204062787 0.324 0.103416435826408],...
    'String',settings,...
    'FontSize',12,...
    'FitBoxToText','off');

end
end

%%
function [Cdata,nCoil,V]= CoilcompressSVD(data,V,tol)
%
%function for coil compression and apply phase for FOV shift
%
%data= COL X COIL X LIN 
%
%Cdata= COLx LINx COIL : compressed data
if(nargin<3)
    tol=0.02;
end
[nCol,nCoil,nLin]=size(data);
D = reshape(permute(data,[1 3 2]),nCol*nLin,nCoil);
if(isempty(V)) 
[~,S,V] = svd(D,'econ');
 nCoil =find(diag(S)/S(1)>tol,1,'last');
 fprintf('Coil compression performed:  %d coils to %d coils with %1.2f tolerance\n',size(V,2),nCoil,tol);
V=V(:,1:nCoil);
else
    nCoil=size(V,2);
end

Cdata = reshape(D*V,nCol,nLin,nCoil);

end
