function [EigenVecs, EigenVals] = kernelEig3D(kernel, imSize,nMaps,crop_val)
% [EigenVecs, EigenVals] = kernelEig3D(kernel, imSize,nMaps,crop_val)
% Function computes the ESPIRiT step II -- eigen-value decomposition of a
% k-space kernel in image space. Kernels should be computed with dat2Kernel
% and then cropped to keep those corresponding to the data space.
%
% INPUTS:
%           kernel - k-space kernels computed with dat2Kernel (5D)
%           imSize - The size of the image to compute maps for [sx,sy,sz]
%           nMaps  - number of ESPIRIT maps
%
% OUTPUTS:
%           EigenVecs - Images representing the Eigenvectors. (sx,sy,sz,Num coils,nMaps)
%           EigenVals - Images representing the EigenValues. (sx,sy,sz,nMaps)
%                       The last are the largest (close to 1)
% praveen.ivp
nCha= size(kernel,4);
nMaps=min(nMaps,nCha);
nElem=nCha*(nCha+1)/2; %number of upper triangular elements
utri_mask_0=triu(ones(nCha))>0;
v=1;
[kernel,v]=getCoilRot((kernel));

kernel_grams= kernel2Gram(kernel,imSize);

if(size(kernel,3)==1),imSize(3)=1; end
KERNEL=sinc_interp3D(kernel_grams,imSize);

KERNEL=reshape(KERNEL,prod(imSize),nElem);


EigenVecs = zeros(prod(imSize), nCha, nMaps,class(KERNEL));
EigenVals = zeros(prod(imSize), nMaps,class(KERNEL));

for n=1:prod(imSize)
    mtx=zeros(nCha,class(KERNEL));
    mtx(utri_mask_0) = squeeze(KERNEL(n,:));
    mtx=triu(mtx,1)'+mtx;
    [C,D,~] = svd(mtx,'econ');

     ph = repmat(exp(-1i*angle(C(1,:))),[size(C,1),1]);
     C = v*(C.*ph);
    D = real(diag(D));
    EigenVals(n,:) = D(1:nMaps);
    EigenVecs(n,:,:) = C(:,1:nMaps);
end
EigenVecs=reshape(EigenVecs,[imSize,nCha,nMaps]);
EigenVals=reshape(EigenVals,[imSize,nMaps]);



%crop?
mask=EigenVals>=crop_val;
EigenVecs=EigenVecs.*reshape(mask,[imSize,1,nMaps]);


end


function kernel_grams= kernel2Gram(kernels,outSize)

sz=size(kernels);

pad_size=sz(1:3);
pad_size(pad_size~=1)=pad_size(pad_size~=1)*2;
% ifft to image space but only 2 times zeropadding:
kernels=zpad(conj(kernels),[pad_size sz(4:end)]);
kernel_im=fftc(fftc(fftc(kernels,1),2),3);
%fix scaling
kernel_im=kernel_im.*sqrt(prod(pad_size)/prod(sz(1:3)));

kernel_im=permute(reshape(kernel_im,[prod(pad_size) sz(4:end)]),[2 3 1]);

% similar to BART compute gram matrices and store only upper triangular
% matrices
nCha=size(kernels,4);
nVoxels=prod(pad_size);

nElem=nCha*(nCha+1)/2; %number of upper triangular elements
utri_mask=triu(ones(nCha))>0;

kernel_grams=zeros(nElem,nVoxels,class(kernels));
for voxel_idx=1:nVoxels

    temp=kernel_im(:,:,voxel_idx)*kernel_im(:,:,voxel_idx)'; % calucalte gram matrices
    kernel_grams(:,voxel_idx)=temp(utri_mask);
end

%reshape, bring physical dims to front and pad to outSize
kernel_grams=reshape(kernel_grams.',[pad_size nElem]);

end


%%
function [kernel,v]=getCoilRot(kernel)
sz=size(kernel);
nCha=size(kernel,4);
kernel=reshape(permute(kernel,[1 2 3 5 4]),[],nCha);

[~,~,v]=svd(kernel'*kernel);

kernel=kernel*v;

kernel=ipermute(reshape(kernel,[sz(1:3) sz(5) sz(4)] ),[1 2 3 5 4]);

end


function test_triu()
% test Mat2TriVec and TriVec2Mat
Q=complex(randn([100 3]),randn([100 3 ]));
QtQ=repmat(Q'*Q,[1  1 2 2 2]);
q_vec=Mat2TriVec(QtQ);
Q_r=TriVec2Mat(q_vec);
assert(sum(QtQ-Q_r,'all')==0) 
end

function out=Mat2TriVec(in)
% out=Mat2TriVec(in);
% convert Hermitian matrices to triu vector

sz=size(in,1:6);
S=sz(1)*(sz(1)+1)/2;
outSize=[S,sz(3:end)];
out=zeros(S,prod(sz(3:end)),class(in));
triu_mask=triu(ones(sz(1)))>0;

in=reshape(in,[],prod(sz(3:end)));

for idx=1:size(in,2)
out(:,idx)=in(triu_mask(:),idx);
end
out=reshape(out,outSize);
end

function out = TriVec2Mat(in)
% out = TriVec2Mat(in)
% convert upper triangular matrix element in dim1 to hermitian or symmetric
% matrices

%S=n*(n+1)/2
getn = @(S) floor((-1 + sqrt(1 + 8 * S)) / 2);
sz=size(in,1:6);
in=reshape(in,sz(1),prod(sz(2:end)));

outSize=[getn(sz(1)),getn(sz(1)) ,sz(2:end)];
assert(mod(outSize(1),1)==0,'vector is of wrong number of elements');
out=zeros([outSize(1:2) outSize(3:end)],class(in));

triu_mask=triu(ones(outSize(1)))>0;

for idx=1:prod(sz(2:end))
    temp=zeros(outSize(1:2));
temp(triu_mask)=(in(:,idx));
temp=temp+triu(temp,1)';
out(:,:,idx)=temp;
end
out=reshape(out,outSize);
end
