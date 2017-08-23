function SSRB_mMR(sino,name)
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

% For reading a file
%sinogramname = strcat('sino_', name, '.raw');
%fid=fopen(sinogramname,'r');

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
%fclose(fid);
clear Sino2D SIN2D;

SINR = double(SSRBSino);

% Gaps Finder and Gap Filling (using inpaint)
MASK = zeros(Nbins,Nproj);
for u=1:Nslices
  MASK = MASK + SINR(:,:,u); 
end  
ngap=0;
nnogap=0;
for ith=1:Nproj
  for ir=1:Nbins
    if (MASK(ir,ith)<0.01)     
      ngap=ngap+1;    
      MASK(ir,ith)=0.0;
    else
      nnogap=nnogap+1; 
      MASK(ir,ith)=1.0;
    end
  end
end  
for u=1:Nslices
  SIN=SINR(:,:,u);
  SIN(MASK==0.)=nan;
  SIN2=inpaint_nans(SIN,2);
  SSRBSino(:,:,u)=SIN2;
end
clear SINR;

% Smooth 3D data with gaussian kernel
SSRBSino = smooth3(SSRBSino,'gaussian',[3 3 3],0.42466);
% sd of 0.42466 correlates to FWHM of 1
%SSRBSino = smooth3(SSRBSino,'gaussian',[5 5 5],1.5);

newName = strcat('sino_SSRB_', name, '.raw');
fid = fopen(newName,'w');
fwrite(fid,SSRBSino,'float32');
fclose(fid);

% % Transform sinogram do STIR projectiondata
% for u=1:Nslices
%   for v=1:Nproj
%     for w=1:Nbins
%       STIR_proj(w,u,v)=SSRBSino(w,v,u);
%     end
%   end
% end
% stirName = strcat('STIR_', filename, '.s');
% fid = fopen(stirName, 'w');
% fwrite(fid,STIR_proj,'float32');
% fclose(fid);
end
