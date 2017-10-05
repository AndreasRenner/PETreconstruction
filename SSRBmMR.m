function [SSRBsino,maxIndex] = SSRBmMR(sino,maxIndexRef)
% Basic Parameters of Siemens Biograph mMR
Nbins   = 344;         % Number of radial bins (NRAD)
Nproj   = 252;         % Number of projections (NANG)
Nseg    = 121;         % Number of segments
Nrings  = 64;          % Number of Detector-Rings
Nslices = 127;         % Number of slices (2*Nrings -1)
% The Sinogram is arranged in 121 Segments; the first Segment consists of
% the 64 direct planes, 2nd and 3rd Segment of 63 planes each for +1 and
% -1 rings and so on; last and second-last Segment consist of 4 planes
% each as the maximum ring difference is 60; this sums up to a total of
% 4084 planes.
Span    = 1;

% SSRB Algorithm: creates a rebinned direct sinogram  
SSRBSino = zeros(Nbins,Nproj,Nslices,'double');
SSRBsino = zeros(Nbins,Nproj,13,'double');
% INCLINATION of each segment
Offset(1)=0;
Offset(2)=(Span+1)/2;
Offset(3)=(Span+1)/2;
for iseg=4:Nseg
  Offset(iseg)=Offset(iseg-2)+Span;
end
% Number of sinograms in each segment
Nsinos=0;
for iseg=1:Nseg
  NsinoSeg(iseg)=Nrings-Offset(iseg);
  Nsinos=Nsinos+NsinoSeg(iseg);
end

Nsino=0;
for iseg=1:Nseg
  for u=1:NsinoSeg(iseg)
    Nsino  = Nsino+1;
    Sino2D = sino(:,:,Nsino);
    % Alternatively read from file
    %Sino2D = fread(fid, [Nbins, Nproj], '*int16');
    k_SSRB = (2*u-1) + Offset(iseg);
    SIN2D  = double(Sino2D); 
    SSRBSino(:,:,k_SSRB)=SSRBSino(:,:,k_SSRB)+SIN2D; 
  end
end
clear Sino2D SIN2D;

if ~maxIndexRef
  % compress sinogram to 13 slices around max-count-slice
  % smooth sinogram
  sino = smooth3(SSRBSino,'gaussian',[5 5 3],1);
  % find slice with highest number of counts
  a = sum(sum(sino));
  maxIndex = find(a==max(a));
  % select 13 slices around max slice
  j = 1;
  for i=(maxIndex-6):(maxIndex+6)
    if i>0
      SSRBsino(:,:,j) = SSRBSino(:,:,i);
    end
    j = j+1;
  end
  fprintf('Current max index is %i\r',maxIndex);
else
  % select 13 slices around max slice
  j = 1;
  for i=(maxIndexRef-6):(maxIndexRef+6)
    if i>0
      SSRBsino(:,:,j) = SSRBSino(:,:,i);
    end
    j = j+1;
  end
end

end
