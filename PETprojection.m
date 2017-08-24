
function PETprojection(filename)
% This program generates the SSRB sinograms corresponding to the segmented
% Object and Transmission Source. 
%
% INPUTS:
% 1 - File 1: Segmented Transmission Source (Tsource, NAC)
% 2 - File 2: Acquired SSRB sinogram
%
% RES: Resolution in XY of the segmented IMAGES 
% NRAD: Number of RADIAL bins (Sinogram)
% NANG: Number of ANGULAR bins (Snogram)
% Nslices: Number os slices of the Image & Sinogram
%
% OUTPUTS:
% Sinogram of just the source (with reduced scatter/random):
%
% QIMP,MAY 2017. 
% Contact info: jacobo.calgonzalez@meduniwien.ac.at
%
% Modified, AUG 2017. Andreas Renner
% -------------------------------------------------------------

% INPUT PARAMETERS
Nparts  = 40;

RES     = 344;
NRAD    = 344;
NANG    = 252;
Nslices = 127;

delta_ang = 180./NANG;

for i=1:Nparts
  % Opening Mask of the PET image and radon projection of this image
  file1 = strcat('mask_',filename,'_',num2str(i),'.raw');
  file2 = strcat('SSRB_cor_',filename,'_',num2str(i),'.raw');
  file3 = strcat('SSRB_seg_',filename,'_',num2str(i),'.raw');
  fid1  = fopen(file1,'r');
  fid2  = fopen(file2,'r');
  fid3  = fopen(file3, 'w');

  %Reading the Mask
  PET2  = zeros(RES,RES,Nslices,'double');
  PET2d = zeros(RES,RES,Nslices,'double');
  for k=1:Nslices
    PET2(:,:,k)  = fread(fid1,[RES,RES], '*float32');
    PET2d(:,:,k) = double(PET2(:,:,k));
    clear PET2;
    for j=1:RES
      for l=1:RES
        if isnan(PET2d(l,j,k));
          PET2d(l,j,k) = 0.0;
        end
      end
    end
  end
  fclose(fid1);

  % ~= ... not equal; does element by element comparisons
  if(RES ~= NRAD);
    PET2dint = PET2d;
    PET2d    = resize(PET2dint,[NRAD,NRAD,Nslices]);
  end

  PET_proj2 = zeros(NRAD,NANG,Nslices,'double');
  for k=1:Nslices
    theta    = 0:delta_ang:179.9;
    [R2,xp2] = radon(PET2d(:,:,k),theta);    
    if(k==1);
      m   = size(xp2,1);
      mmm = (m-NRAD)/2-0.5;
    end
    for l=1:RES 
      PET_proj2(l,:,k) = R2(l+mmm,:);      
    end
  end

  for k=1:Nslices
    proj2d = PET_proj2(:,:,k);
  end

  factor2 = zeros(NRAD,NANG,Nslices,'double');
  for k=1:Nslices
    suma = PET_proj2(:,:,k);
    factor2(:,:,k) = PET_proj2(:,:,k)./suma;
  end

  for k=1:Nslices
    input(:,:,k)  = fread(fid2, [NRAD,NANG], '*float32');
    output(:,:,k) = input(:,:,k).*factor2(:,:,k);
    for j=1:NRAD
      for l=1:NANG
        if isnan(output(j,l,k));
          output(j,l,k) = 0.0;
        end
      end
    end
  end
  fwrite(fid3,output,'float32');
  fclose(fid2);
  fclose(fid3);
end

end
