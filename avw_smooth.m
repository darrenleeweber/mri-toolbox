function avw = avw_smooth(avw,fwhm)

% avw_smooth - Guassian smoothing
%
% Usage: avw = avw_smooth(avw,fwhm)
%
% avw is the Analyze struct returned by avw_read
% fwhm is an integer indicating the voxels to convolve with the Gaussian kernel
% (default = 3, a 3x3x3 convolution)
%

if ~exist('fwhm','var'), fwhm = 3; end
if isempty(fwhm), fwhm = 3; end

fprintf('...gaussian smoothing...'); tic;

avw.img = smooth3(avw.img,'gaussian',fwhm);

t = toc; fprintf('done (%5.2f sec)\n',t);

return