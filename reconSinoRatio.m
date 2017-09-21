function reconSinoRatio(nameBlank,nameTrans,reconMethod)
% Reconstruct ratio between Blanc- and Transmission-Scan

filenameBlank = strcat('SSRB_complete_',  nameBlank,'.raw');
filenameTrans = strcat('SSRB_complete_',  nameTrans,'.raw');
filenameMask  = strcat('mask_SSRB_ratio_',nameTrans,'.raw');

% Basic Parameters of Siemens Biograph mMR
Nbins   = 344;         % Number of radial bins (NRAD)
Nproj   = 252;         % Number of projections (NANG)
Nslices = 127;         % Number of slices (2*Nrings -1)
% SSRB Algorithm: creates a rebinned direct sinogram  
% SSRBSino = zeros(Nbins,Nproj,Nslices,'double');
%   -> size(SSRB,1) = 344
%   -> size(SSRB,2) = 252
%   -> size(SSRB,3) = 127

%load('SSRB_Blanc.mat');
%load('SSRB_Trans.mat');
SinoBlank = zeros(Nbins,Nproj,Nslices,'double');
SinoTrans = zeros(Nbins,Nproj,Nslices,'double');
SinoMask  = zeros(Nbins,Nproj,Nslices,'double');
fid1      = fopen(filenameBlank,'r');
fid2      = fopen(filenameTrans,'r');
fid3      = fopen(filenameMask, 'r');
for i=1:Nslices
  Blank2D = fread(fid1,[Nbins,Nproj],'float32');
  Trans2D = fread(fid2,[Nbins,Nproj],'float32');
  Mask2D  = fread(fid3,[Nbins,Nproj],'float32');
  SinoBlank(:,:,i) = SinoBlank(:,:,i) + Blank2D;
  SinoTrans(:,:,i) = SinoTrans(:,:,i) + Trans2D;
  SinoMask(:,:,i)  = SinoMask(:,:,i)  + Mask2D;
end
fclose(fid1);
fclose(fid2);
fclose(fid3);

% KernelSize should be 2*ceil(2*sigma)+1
%SinoBlank = smooth3(SinoBlank,'gaussian',[9 9 7],2);
%SinoTrans = smooth3(SinoTrans,'gaussian',[9 9 7],2);
%SinoTrans = smooth3(SinoTrans,'gaussian',[5 5 5],1);
%SinoBlank = smooth3(SinoBlank,'gaussian',[5 5 5],1);
SinoTrans = smooth3(SinoTrans,'gaussian',[3 3 3],0.42466);
SinoBlank = smooth3(SinoBlank,'gaussian',[3 3 3],0.42466);

SinoRatio = zeros(Nbins,Nproj,Nslices,'double');
%SinoRatio = log(SinoBlank)-log(SinoTrans);

for u=1:Nbins
  for v=1:Nproj
    for w=1:Nslices
      if SinoTrans(u,v,w)<=0 || SinoBlank(u,v,w)<=0
        SinoRatio(u,v,w) = 0.0;
      elseif SinoTrans(u,v,w)<0.1 % SSRB_Blank(u,v,w)<4 &&
        % Factor 5 -> pixel values correspond to attenuation per cm
        % -> is done in final reconstruction
        SinoRatio(u,v,w) = log(SinoBlank(u,v,w));
      else
        % Factor 5 -> pixel values correspond to attenuation per cm
        % -> is done in final reconstruction
        SinoRatio(u,v,w) = log(SinoBlank(u,v,w)/(SinoTrans(u,v,w)));
      end
    end
  end
end

%SinoRatio=smooth3(SinoRatio,'gaussian',[9 9 7],2);
SinoRatio=smooth3(SinoRatio,'gaussian',[3 3 3],0.42466);

name = strcat('SSRB_333ratio333_',nameTrans,'.raw');
fid = fopen(name,'w');
fwrite(fid,SinoRatio,'float32');
fclose(fid);

% ****************************
% ToDo:
% multiplication with mask
% ****************************

if strcmp(reconMethod, 'FBP')
  theta = 180./Nproj;
  for i=1:Nslices
    recon(:,:,i) = 4.83559*iradon(SinoRatio(:,:,i),theta,Nbins);
  end
  name = strcat('recon_333ratio333_',nameTrans,'.raw');
  fid  = fopen(name, 'w');
  fwrite(fid,recon,'float32');
  fclose(fid);
elseif strcmp(reconMethod, 'OSEM')
  recon = 4.83559*OSEM_Recon(SinoRatio,nameTrans);
  name = strcat('OSEM_333ratio333_',nameTrans,'.raw');
  fid  = fopen(name, 'w');
  fwrite(fid,recon,'float32');
  fclose(fid);
else
  fprintf('Wrong reconMethod - options are: "FBP" and "OSEM"!');
end

end