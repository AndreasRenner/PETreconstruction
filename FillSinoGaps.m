function Sino = FillSinoGaps(filename,maskname)
% Function to inpaint_nans in a ROI defined by a mask

% Basic Parameters of Siemens Biograph mMR
Nbins   = 344;         % Number of radial bins (NRAD)
Nproj   = 252;         % Number of projections (NANG)
Nslices = 127;         % Number of slices (2*Nrings -1)

Sino = zeros(Nbins,Nproj,Nslices,'double');
Mask = zeros(Nbins,Nproj,Nslices,'double');
%recon= zeros(Nbins,Nbins,Nslices,'double');

SE = strel('disk',30);

fid1 = fopen(filename,'r');
fid2 = fopen(maskname,'r');
for i=1:Nslices
  Sino2D = fread(fid1,[Nbins,Nproj],'float32');
  Mask2D = fread(fid2,[Nbins,Nproj],'float32');
  Mask2D = imdilate(Mask2D,SE);
  Mask2D = imgaussfilt(Mask2D,10);
  for j=1:Nproj
    for k=1:Nbins
      if Mask2D(k,j) && Sino2D(k,j)<=0
        Sino2D(k,j) = nan;
      end
    end
  end
  Sino2D = inpaint_nans(Sino2D,2);
  Sino(:,:,i) = Sino(:,:,i) + Sino2D;
  if i>10 % to avoide noise in top of image
    Mask(:,:,i) = Mask(:,:,i) + Mask2D;
  end
end
fclose(fid1);
fclose(fid2);

% Mask Sinogram and write result to file
Sino = Sino.*Mask;
%Sino = smooth3(Sino,'gaussian',[9 9 7],2);
%name = strcat('gapfilled<1_ratio_02TransPhantom2_31-08.raw');
%fid = fopen(name, 'w');
%fwrite(fid,Sino,'float32');
%fclose(fid);

% Reconstruct Sinogram
%theta = 180./Nproj;
%for i=1:Nslices
%  recon(:,:,i)=iradon(Sino(:,:,i),theta,Nbins);
%end
% Write reconstruction to file
%name = strcat('recon_gapfilled<1_ratio_02TransPhantom2_31-08.raw');
%fid = fopen(name, 'w');
%fwrite(fid,recon,'float32');
%fclose(fid);
end