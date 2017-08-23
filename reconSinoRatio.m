function reconSinoRatio(SSRB_Blank,SSRB_Trans,filename)
% Reconstruct ratio between Blanc- and Transmission-Scan

% Basic Parameters of Siemens Biograph mMR
% Nbins   = 344;         % Number of radial bins (NRAD)
% Nproj   = 252;         % Number of projections (NANG)
% Nslices = 127;         % Number of slices (2*Nrings -1)
% SSRB Algorithm: creates a rebinned direct sinogram  
% SSRBSino = zeros(Nbins,Nproj,Nslices,'double');
%   -> size(SSRB,1) = 344
%   -> size(SSRB,2) = 252
%   -> size(SSRB,3) = 127

%load('SSRB_Blanc.mat');
%load('SSRB_Trans.mat');

% KernelSize should be 2*ceil(2*sigma)+1
%SSRB_Blank=smooth3(SSRB_Blank,'gaussian',[9 9 7],2);
%SSRB_Trans=smooth3(SSRB_Trans,'gaussian',[9 9 7],2);
%SSRB_Trans=smooth3(SSRB_Trans,'gaussian',[5 5 5],1);
%SSRB_Blank=smooth3(SSRB_Blank,'gaussian',[5 5 5],1);
SSRB_Trans=smooth3(SSRB_Trans,'gaussian',[3 3 3],0.42466);
SSRB_Blank=smooth3(SSRB_Blank,'gaussian',[3 3 3],0.42466);

SSRB_Ratio=zeros(size(SSRB_Blank,1),size(SSRB_Blank,2),size(SSRB_Blank,3),'double');
%SSRB_Ratio=5*(log(SSRB_Blank)-log(SSRB_Trans));

for u=1:size(SSRB_Blank,1)
  for v=1:size(SSRB_Blank,2)
    for w=1:size(SSRB_Blank,3)
      if SSRB_Trans(u,v,w)<=0 || SSRB_Blank(u,v,w)<=0
        SSRB_Ratio(u,v,w)=0.0;
      elseif SSRB_Trans(u,v,w)<0.1 % SSRB_Blank(u,v,w)<4 &&
        % Factor 5 -> pixel values correspond to attenuation per cm
        SSRB_Ratio(u,v,w)=5*log(SSRB_Blank(u,v,w));
      else
        % Factor 5 -> pixel values correspond to attenuation per cm
        SSRB_Ratio(u,v,w)=5*log(SSRB_Blank(u,v,w)./(SSRB_Trans(u,v,w)));
      end
    end
  end
end

%SSRB_Ratio=smooth3(SSRB_Ratio,'gaussian',[9 9 7],2);
SSRB_Ratio=smooth3(SSRB_Ratio,'gaussian',[3 3 3],0.42466);

name = strcat('SSRB_333ratio_',filename,'.raw');
fid = fopen(name,'w');
fwrite(fid,SSRB_Ratio,'float32');
fclose(fid);

theta = 180./size(SSRB_Ratio,2);
for i=1:size(SSRB_Blank,3)
  recon(:,:,i)=iradon(SSRB_Ratio(:,:,i),theta,size(SSRB_Ratio,1));
end

name = strcat('recon_333ratio_',filename,'.raw');
fid = fopen(name, 'w');
fwrite(fid,recon,'float32');
fclose(fid);
end