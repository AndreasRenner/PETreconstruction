function segmentTube(blankname,transname,VoxelPerSegment)
% INPUTS:               OUTPUTS:
% -> sino_cor_blank     -> mask and SSRB_mask
% -> sino_cor_trans     -> sino_seg_blank and trans
% -> recon_blank        -> complete_blank and trans + totalmask

% Basic Parameters of Siemens Biograph mMR
Nbins    = 344;         % Number of radial bins (NRAD)
Nproj    = 252;         % Number of projections (NANG)
Nslices  = 127;         % Number of slices (2*Nrings -1)
deltaAng = 180./Nproj;

Nparts   = 160;         % Number of sinogram parts within one scan

% Read Segment-Mask
segMask1 = zeros(Nbins,Nbins,Nslices);
segMask2 = zeros(Nbins,Nbins,Nslices);
fid1 = fopen('Mask1.raw','r');
fid2 = fopen('Mask2.raw','r');
for i=1:Nslices
  mask1 = fread(fid1,[Nbins,Nbins],'float32');
  segMask1(:,:,i) = segMask1(:,:,i) + mask1;
  mask2 = fread(fid2,[Nbins,Nbins],'float32');
  segMask2(:,:,i) = segMask2(:,:,i) + mask2;
end
fclose(fid1);
fclose(fid2);

% Initialaize Sum of all Mask-Projections
sinoMaskTotal = zeros(Nbins,Nproj,Nslices,'double');

for i=1:Nparts
  % Read corrected sinograms and reconstruction of tube segments
  SinoBlank = zeros(Nbins,Nproj,Nslices,'double');
  SinoTrans = zeros(Nbins,Nproj,Nslices,'double');
  %sinoMask  = zeros(Nbins,Nproj,Nslices,'double');
  name1 = strcat('SSRB_cor_',blankname,'_',num2str(i),'.raw');
  name2 = strcat('SSRB_cor_',transname,'_',num2str(i),'.raw');
  %name3 = strcat('SSRB_mask_600_',blankname,'_',num2str(i),'.raw');
  fid1 = fopen(name1,'r');
  fid2 = fopen(name2,'r');
  %fid3 = fopen(name3,'r');
  for l=1:Nslices
    Blank2D = fread(fid1,[Nbins,Nproj],'float32');
    SinoBlank(:,:,l) = SinoBlank(:,:,l) + Blank2D;
    Trans2D = fread(fid2,[Nbins,Nproj],'float32');
    SinoTrans(:,:,l) = SinoTrans(:,:,l) + Trans2D;
    %Recon2D = fread(fid3,[Nbins,Nproj],'float32');
    %sinoMask(:,:,l) = sinoMask(:,:,l) + Recon2D;
  end
  fclose(fid1);
  fclose(fid2);
  %fclose(fid3);
  
  % Reconstruction of tube segments using FBP
  SinoTotal = SinoBlank+SinoTrans;
  for u=1:Nslices
    ReconSegm(:,:,u) = iradon(SinoTotal(:,:,u),deltaAng,Nbins);
  end
  %ReconSegm = smooth3(ReconSegm,'gaussian',[3 3 3],0.42466);
    
  % Segmentation of the reconstructed tube segments
  mask    = zeros(Nbins,Nbins,Nslices,'double');
  tmpmask = zeros(Nbins,Nbins,Nslices,'double');
  tmpmask = ReconSegm;
  flatmask = reshape(tmpmask,Nbins*Nbins*Nslices,1);
  [sortmask,index] = sort(flatmask,'descend');
  clear sortmask;
  flatmask = zeros(Nbins*Nbins*Nslices,1);
  for u=1:VoxelPerSegment
    flatmask(index(u))=1;
  end
  tmpmask = reshape(flatmask,[Nbins,Nbins,Nslices]);
  mask = tmpmask;
  
  %mask = smooth3(mask,'gaussian',[3 3 3],0.42466);
  
  % Apply segment-mask to segmentation of reconstructed tube segments
  if mod(i,2)
    mask = mask.*segMask1;
  else
    mask = mask.*segMask2;
  end

  % Write segmentation result (-> mask) to file
  name = strcat('mask_',num2str(VoxelPerSegment),'_',blankname,'_',num2str(i),'.raw');
  fid = fopen(name, 'w');
  fwrite(fid,mask,'float32');
  fclose(fid);
  
  % Reproject mask to get it in sinogram-space
  PET_proj2 = zeros(Nbins,Nproj,Nslices,'double');
  for k=1:Nslices
    theta    = 0:deltaAng:179.9;
    [R2,xp2] = radon(mask(:,:,k),theta);    
    if(k==1);
      m   = size(xp2,1);
      mmm = (m-Nbins)/2-0.5;
    end
    for l=1:Nbins 
      PET_proj2(l,:,k) = R2(l+mmm,:);      
    end
  end  
  sinoMask = zeros(Nbins,Nproj,Nslices,'double');
  for k=1:Nslices
    suma = PET_proj2(:,:,k);
    sinoMask(:,:,k) = PET_proj2(:,:,k)./suma;
    for j=1:Nbins
      for l=1:Nproj
        if isnan(sinoMask(j,l,k))
          sinoMask(j,l,k)=0.0;
        end
      end
    end
  end
  
  sinoMaskTotal = sinoMaskTotal + sinoMask;
  
  % *******************************************
  % ToDo:
  % add summation over SSRB_cor to create reference for scatter
  % correction -> how much does it improve?
  % *******************************************
  
  % Write sino mask to file
  name = strcat('SSRB_mask_',num2str(VoxelPerSegment),'_',blankname,'_',num2str(i),'.raw');
  fid = fopen(name, 'w');
  fwrite(fid,sinoMask,'float32');
  fclose(fid);
end
  
% Write Total Sino-Mask to correct for overlap of sinograms
sinoMaskTotal = sinoMaskTotal + 1;
name = strcat('SSRB_mask_TOTAL_',num2str(VoxelPerSegment),'_',blankname,'.raw');
fid = fopen(name,'w');
fwrite(fid,sinoMaskTotal,'float32');
fclose(fid);

for i=1:Nparts
  % Read corrected sinograms and reconstruction of tube segments
  SinoBlank = zeros(Nbins,Nproj,Nslices,'double');
  SinoTrans = zeros(Nbins,Nproj,Nslices,'double');
  sinoMask  = zeros(Nbins,Nproj,Nslices,'double');
  name1 = strcat('SSRB_cor_',blankname,'_',num2str(i),'.raw');
  name2 = strcat('SSRB_cor_',transname,'_',num2str(i),'.raw');
  name3 = strcat('SSRB_mask_',num2str(VoxelPerSegment),'_',blankname,'_',num2str(i),'.raw');
  fid1 = fopen(name1,'r');
  fid2 = fopen(name2,'r');
  fid3 = fopen(name3,'r');
  for l=1:Nslices
    Blank2D = fread(fid1,[Nbins,Nproj],'float32');
    SinoBlank(:,:,l) = SinoBlank(:,:,l) + Blank2D;
    Trans2D = fread(fid2,[Nbins,Nproj],'float32');
    SinoTrans(:,:,l) = SinoTrans(:,:,l) + Trans2D;
    Mask2D = fread(fid3,[Nbins,Nproj],'float32');
    sinoMask(:,:,l) = sinoMask(:,:,l) + Mask2D;
  end
  fclose(fid1);
  fclose(fid2);
  fclose(fid3);
  
  % Apply reprojected mask to blank- and trans-sinogram
  SinoBlank = SinoBlank.*sinoMask./sinoMaskTotal;
  SinoTrans = SinoTrans.*sinoMask./sinoMaskTotal;
  name1 = strcat('SSRB_seg_',num2str(VoxelPerSegment),'_',blankname,'_',num2str(i),'.raw');
  name2 = strcat('SSRB_seg_',num2str(VoxelPerSegment),'_',transname,'_',num2str(i),'.raw');
  fid1 = fopen(name1,'w');
  fid2 = fopen(name2,'w');
  fwrite(fid1,SinoBlank,'float32');
  fwrite(fid2,SinoTrans,'float32');
  fclose(fid1);
  fclose(fid2);
end

% Add up all parts to get Sinograms of whole acquisition
SinoBlank = zeros(Nbins,Nproj,Nslices,'double');
SinoTrans = zeros(Nbins,Nproj,Nslices,'double');
TotalMask = zeros(Nbins,Nbins,Nslices,'double');
for i=1:Nparts
  name1 = strcat('SSRB_seg_',num2str(VoxelPerSegment),'_',blankname,'_',num2str(i),'.raw');
  name2 = strcat('SSRB_seg_',num2str(VoxelPerSegment),'_',transname,'_',num2str(i),'.raw');
  name3 = strcat('mask_',num2str(VoxelPerSegment),'_',blankname,'_',num2str(i),'.raw');
  fid1  = fopen(name1,'r');
  fid2  = fopen(name2,'r');
  fid3  = fopen(name3,'r');
  for l=1:Nslices
    Blank2D = fread(fid1,[Nbins,Nproj],'float32');
    Trans2D = fread(fid2,[Nbins,Nproj],'float32');
    Mask2D  = fread(fid3,[Nbins,Nbins],'float32');
    for j=1:Nbins
      for k=1:Nproj
        if isnan(Blank2D(j,k));
          Blank2D(j,k) = 0.0;
        end
        if isnan(Trans2D(j,k));
          Trans2D(j,k) = 0.0;
        end
      end
    end
    SinoBlank(:,:,l) = SinoBlank(:,:,l) + Blank2D;
    SinoTrans(:,:,l) = SinoTrans(:,:,l) + Trans2D;
    TotalMask(:,:,l) = TotalMask(:,:,l) + Mask2D;
  end
  fclose(fid1);
  fclose(fid2);
  fclose(fid3);
end

% Write Sinograms to file
name1 = strcat('SSRB_complete_',num2str(VoxelPerSegment),'_',blankname,'.raw');
name2 = strcat('SSRB_complete_',num2str(VoxelPerSegment),'_',transname,'.raw');
name3 = strcat('totalMask_',num2str(VoxelPerSegment),'_',blankname,'.raw');
fid1  = fopen(name1,'w');
fid2  = fopen(name2,'w');
fid3  = fopen(name3,'w');
fwrite(fid1,SinoBlank,'float32');
fwrite(fid2,SinoTrans,'float32');
fwrite(fid3,TotalMask,'float32');
fclose(fid1);
fclose(fid2);
fclose(fid3);

totalCounts = sum(sum(sum(SinoTrans)));
fprintf('Total Counts in Trans-Sinogram: %u\r',totalCounts);

reconblank = strcat(num2str(VoxelPerSegment),'_',blankname);
recontrans = strcat(num2str(VoxelPerSegment),'_',transname);
reconSinoRatio(reconblank,recontrans,'FBP');
reconSinoRatio(reconblank,recontrans,'OSEM');

end