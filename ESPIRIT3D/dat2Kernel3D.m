function [kernels, S] = dat2Kernel3D(data, kSize,threshold)
% kernel = dat2Kernel(data, kSize,thresh)
%
% Function to perform k-space calibration step for ESPIRiT and create
% k-space kernels. Only works for 2D multi-coil images for now.  
% 
% Inputs: 
%       data - calibration data [kx,ky,kz,coils]
%       kSize - size of kernel (for example kSize=[6,6,6])
%       threshold - smallest normalized singular value to keep (0,1]
%
% Outputs: 
%       kernel - k-space kernels matrix (not cropped), which correspond to
%                the basis vectors of overlapping blocks in k-space
%       S      - (Optional parameter) The singular vectors of the
%                 calibration matrix
%

[sx,sy,sz,nc] = size(data);

kSize = min(kSize,[sx,sy,sz]);


tmp = im2row3D(data,kSize); [tsx,tsy,tsz] = size(tmp);
A = reshape(tmp,tsx,tsy*tsz);
 A=A'*A;

[~,S,V] = svd(A,'econ');  
S = sqrt(abs(diag(S)));

% pick kernels now 
figure,subplot(211),plot(S),hold on, stem(sum(S>S(1)*threshold),S(1)); %debug plot
subplot(212),imagesc(abs(V)),title('singualr vectors'),colormap('gray')
S=real(S);
S=S(S>(S(1)*threshold));
fprintf('picking %d kernels with threshold %.4f of largest singular value %d\n',length(S),threshold,S(1));

V=V(:,1:length(S));
kernels = reshape(V,kSize(1),kSize(2),kSize(3),nc,size(V,2));
end


function ImRow=im2row3D(im,winSize)

%INPUTS:
% im: 4D matrix [CHAxLINxPARxNCHA]
%KernelSize: 3 element vector
%OUTPUTS:
%imRow:[Nkernels x KernelPoints x NCHA]



[sx,sy,sz,NCha] = size(im);
Nkernel=(sx-winSize(1)+1)*(sy-winSize(2)+1)*(sz-winSize(3)+1);
ImRow = zeros((sx-winSize(1)+1)*(sy-winSize(2)+1)*(sz-winSize(3)+1),prod(winSize),NCha,class(im));
count=0;
for z=1:winSize(3)
for y=1:winSize(2)
    for x=1:winSize(1)
        count = count+1;
        ImRow(:,count,:) = reshape(im(x:sx-winSize(1)+x,y:sy-winSize(2)+y,z:sz-winSize(3)+z,:),...
            Nkernel,1,NCha);
    end
end
end

end
