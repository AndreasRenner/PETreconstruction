function createRefSino(filename)
% Reconstruction of 'filename'

% Basic Parameters of Siemens Biograph mMR
Nbins   = 344;         % Number of radial bins (NRAD)
Nproj   = 252;         % Number of projections (NANG)
Nslices = 127;         % Number of slices (2*Nrings -1)

delta_ang = 180/Nproj;

% Allocate memory for masks
maskWater  = zeros(Nbins,Nbins,Nslices,'double');
maskAir    = zeros(Nbins,Nbins,Nslices,'double');
maskTeflon = zeros(Nbins,Nbins,Nslices,'double');

% load all Masks
nameMaskWater  = strcat(filename,'_Mask(>0.03).raw');
nameMaskAir    = strcat(filename,'_AirMask(0.022-0.075).raw'); 
nameMaskTeflon = strcat(filename,'_TeflonMask(>0.12).raw');
nameMaskTotal  = strcat(filename,'_MaskTotal.raw');
fid1  = fopen(nameMaskWater,'r');
fid2  = fopen(nameMaskAir,'r');
fid3  = fopen(nameMaskTeflon,'r');
for i=1:Nslices
  MaskWater2D  = fread(fid1,[Nbins,Nbins],'float32');
  MaskAir2D    = fread(fid2,[Nbins,Nbins],'float32');
  MaskTeflon2D = fread(fid3,[Nbins,Nbins],'float32');
  maskWater(:,:,i)  = maskWater(:,:,i)  + MaskWater2D;
  maskAir(:,:,i)    = maskAir(:,:,i)    + MaskAir2D;
  maskTeflon(:,:,i) = maskTeflon(:,:,i) + MaskTeflon2D;
end
fclose(fid1);
fclose(fid2);
fclose(fid3);

% convert mask to binary format
maskWater  = logical(maskWater);
maskAir    = logical(maskAir);
maskTeflon = logical(maskTeflon);

%% Perform post-processing of mask
% create spherical structure element with radius 10
SE = strel('disk',10);
% close the main mask (maskWater)
maskWater = imclose(maskWater,SE);
% smooth boundaries of maskWater using erode and dilate
maskWater = imerode(maskWater,SE);
maskWater = imdilate(maskWater,SE);
% save processed mask for visual check
name = strcat('new_',nameMaskWater);
fid = fopen(name, 'w');
fwrite(fid,maskWater,'float32');
fclose(fid);
%
SEs = strel('disk',4);
maskTeflon = imerode(maskTeflon,SEs);
maskTeflon = imdilate(maskTeflon,SEs);
% save processed mask for visual check
name = strcat('new_',nameMaskTeflon);
fid = fopen(name, 'w');
fwrite(fid,maskTeflon,'float32');
fclose(fid);
%
SEs = strel('disk',4);
SEb = strel('disk',5);
maskAir = imerode(maskAir,SEb);
maskAir = imdilate(maskAir,SEs);
% save processed mask for visual check
name = strcat('new_',nameMaskAir);
fid = fopen(name, 'w');
fwrite(fid,maskAir,'float32');
fclose(fid);

maskTotal = maskWater * 0.096 / 4.83559;
for i=1:Nbins
  for j=1:Nbins
    for k=1:Nslices
      if maskTeflon(i,j,k)
        maskTotal(i,j,k) = 0.184 / 4.83559; % value according to NIST database
      elseif maskAir(i,j,k)
        maskTotal(i,j,k) = 0.0;
      end      
    end
  end
end
% save processed mask for visual check
name = strcat('new_',nameMaskTotal);
fid = fopen(name, 'w');
fwrite(fid,maskTotal,'float32');
fclose(fid);

%% Reproject maskTotal
maskReproj = zeros(Nbins,Nproj,Nslices,'double');
for k=1:Nslices
  theta    = 0:delta_ang:179.9;
  [R2,xp2] = radon(maskTotal(:,:,k),theta);    
  if(k==1);
    m   = size(xp2,1);
    mmm = (m-Nbins)/2-0.5;
  end
  for l=1:Nbins 
    maskReproj(l,:,k) = R2(l+mmm,:);      
  end
end
% save reprojected mask for visual check
name = strcat('SSRB_',nameMaskTotal);
fid = fopen(name, 'w');
fwrite(fid,maskReproj,'float32');
fclose(fid);

end