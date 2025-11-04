function out=sinc_interp3D(data,outsize)
% out=sinc_interp3D(data,outsize)
% test:  as(sinc_interp3D(phantom3d(64),[1,1,1]*256))
assert(length(outsize)==3,'output dimension has to have 3 elements')
out=fftc(fftc(fftc(data,1),2),3);
pad_size=[outsize,size(data,4:max(5,ndims(data)))];
out=zpad(out,pad_size);
out=ifftc(ifftc(ifftc(out,1),2),3);

%fix dumb scaling
out=out.*sqrt(prod(outsize(1:3))/prod(size(data,1:3)));

end

