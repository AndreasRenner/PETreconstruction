function SSRB_add_Scans(filename)

% Basic Parameters of Siemens Biograph mMR
Nbins   = 344;         % Number of radial bins (NRAD)
Nproj   = 252;         % Number of projections (NANG)
Nslices = 127;         % Number of slices (2*Nrings -1)

% Number of sinogram parts within one scan
Nparts = 160;
parts  = Nparts/20;

% Add time from filling of Pellet to start of total scan
if     strcmp(filename, '05BlankFast')
  decayCor = 3471;
  Nscans = 10;
  cd '05BlankFast'
elseif strcmp(filename, '02TransPhantom1')
  decayCor = 2234;
  Nscans = 2;
  cd '02TransPhantom1'
elseif strcmp(filename, '03TransPhantom2')
  decayCor = 2630;
  Nscans = 2;
elseif strcmp(filename, '04TransRQ')
  decayCor = 3039;
  Nscans = 2;
elseif strcmp(filename, '01Blank')
  decayCor = 930;
  Nscans = 15;
  dataFolder = '/media/andreas/CorneaOCT_StudyData2/01Blank_0909';
  mainFolder = cd(dataFolder);
  nameSegMask1 = 'Mask1.raw';
  nameSegMask2 = 'Mask2.raw';
elseif strcmp(filename, '02TransPhantom2')
  decayCor = 2910;
  Nscans = 6;
  %dataFolder = '/media/andreas/CorneaOCT_StudyData2/Backup/02TransPhantom2';
  %mainFolder = cd(dataFolder);
elseif strcmp(filename, '03TransPhantom1')
  decayCor = 3668;
  Nscans = 7;
end

% Load all Files into a list
name = strcat('sino_SSRB_',filename,'_*_*.raw');
d = dir(name);
filenames = {d.name};

%cd(mainFolder);

difference = zeros(Nscans*20,1);
k=0;
for i=1:Nscans
  name = strcat('minlistScan',num2str(i),filename);
  load(name);
  for j=1:20
    k=k+1;
    difference(k) = minlist(j+1)-minlist(j);
  end
end

halflifeFDG= 6586.2;   % halflife of FDG in [s]
decayF     = zeros(Nparts*Nscans,1);
timeF      = zeros(Nparts*Nscans,1);
timespan   = zeros(Nparts*Nscans,1);
totalCorF  = zeros(Nparts*Nscans,1);
partnumber = zeros(Nparts*Nscans,1);

if strcmp(filename, '01Blank')
  segMask1   = zeros(Nbins,Nbins,Nslices);
  segMask2   = zeros(Nbins,Nbins,Nslices);
  % Read Segment-Mask
  fid1 = fopen(nameSegMask1,'r');
  fid2 = fopen(nameSegMask2,'r');
  for i=1:Nslices
    mask1 = fread(fid1,[Nbins,Nbins],'float32');
    segMask1(:,:,i) = segMask1(:,:,i) + mask1;
    mask2 = fread(fid2,[Nbins,Nbins],'float32');
    segMask2(:,:,i) = segMask2(:,:,i) + mask2;
  end
  fclose(fid1);
  fclose(fid2);
end

% (1) DECAY CORRECTION
% calculate decay correction factor for each part
for i=1:(Nparts*Nscans)
  tmpfile = cell2mat(filenames(i));
  fileinfo = strsplit(tmpfile,'_');
  numinfo = cell2mat(fileinfo(5));
  number = strsplit(numinfo,'.');
  partnumber(i) = str2num(cell2mat(number(1)));
  timestamp = str2num(cell2mat(fileinfo(4)))/1000;
  decayF(i) = 2^(-(decayCor+timestamp)/halflifeFDG);
end

% use the sorted counter as reference for the time correction
% obtained from difference(i)
[sortdecay,index] = sort(decayF,'descend');
clear sortdecay;

% (2) SCAN-TIME CORRECTION
% calculate scan-time correction factor for each part
for i=1:(Nscans*20)
  for j=1:parts
    timespan(parts*i-(j-1)) = difference(i)/parts;
  end
end
for i=1:(Nparts*Nscans)
  meantime = 1000;
  timeF(i) = meantime/timespan(index(i));
end

% (3) SCAN-NUMBER CORRECTION
scanNumCor = 10/Nscans;

% calculate corrected sinograms
theta = 180./Nproj;
for i=1:Nparts
  SSRBSino = zeros(Nbins,Nproj,Nslices,'double');
  
  % Find indices of pairs
  indices  = find(partnumber==i); 
  
  %cd(dataFolder);
  
  for j=1:Nscans
    k = indices(j);
    totalCorF(k) = scanNumCor*timeF(k)/decayF(k);
    %fprintf('Total Cor. Faktor for part %u, %u is: %u\r',i,j,totalCorF(k));
    name = cell2mat(filenames(k));
    fid  = fopen(name,'r');
    for l=1:Nslices
      Sino2D = fread(fid,[Nbins,Nproj],'float32');
      SIN2D  = double(Sino2D)*totalCorF(k);
      SSRBSino(:,:,l) = SSRBSino(:,:,l) + SIN2D;
    end
    fclose(fid);
  end
  
  %cd(mainFolder);
  
  % Fil Crystal Gaps (using inpaint)
  SINR = double(SSRBSino);
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
  SSRBSino = smooth3(SSRBSino,'gaussian',[3 3 3],0.42466);
  
%  D E P R I C A T E D
%  -> destroys quantitation
%  % (4) SENSITIVITY CORRECTION
%  % calculate total number of cor. counts
%  totalcounts = sum(sum(sum(SSRBSino)));
%  % apply correction factor
%  sensitivityCor = 200000000./totalcounts;
%  %fprintf('Total Counts for part %u is: %u; corF: %u\r',i,totalcounts,sensitivityCor);
%  SSRBSino = SSRBSino*sensitivityCor;
  
  % Write fully corrected sinogram to file
  name = strcat('SSRB_cor_',filename,'_',num2str(i),'.raw');
  fid = fopen(name,'w');
  fwrite(fid,SSRBSino,'float32');
  fclose(fid);
  
  if strcmp(filename, '01Blank')
    % RECONSTRUCTION
    % Reconstruction of Part using FBP
    recon = zeros(Nbins,Nbins,Nslices,'double');
    for u=1:Nslices
      recon(:,:,u) = iradon(SSRBSino(:,:,u),theta,Nbins);
    end
    % Apply Segment-Mask to recon to get final result
    if mod(i,2)
      recon = recon.*segMask1;
    else
      recon = recon.*segMask2;
    end
    % Write reconstructed image to file
    name = strcat('recon_', filename,'_',num2str(i),'.raw');
    fid = fopen(name, 'w');
    fwrite(fid,recon,'float32');
    fclose(fid);
  
    % SEGMENTATION
    % Segmentation of the reconstructed part
    mask = zeros(Nbins,Nbins,Nslices,'double');
    tmpmask = zeros(Nbins,Nbins,Nslices,'double');
    %tmpmask = zeros(Nbins,Nbins,49,'double');
    %counter = 1;
    %step = floor(2*(41-i));
    %for u = (step-1):(step+47)
      % Here one could add additional segmentation regarding radius
      % to distinguish between emission and transmission
    %  tmpmask(:,:,counter) = recon(:,:,u);
    %  counter = counter+1;
    %end
    tmpmask = recon;
    flatmask = reshape(tmpmask,Nbins*Nbins*Nslices,1);
    %flatmask = reshape(tmpmask,Nbins*Nbins*49,1);
    [sortmask,index] = sort(flatmask,'descend');
    clear sortmask;
    flatmask = zeros(Nbins*Nbins*Nslices,1);
    %flatmask = zeros(Nbins*Nbins*49,1);
    % Find 7200 voxel/round with biggest intensity in tmpmask
    % -> corresponds to diameter of 10 mm (=Pellet+2 mm in every
    %    direction -> max range of positrons)
    for u=1:(7200/parts)
      flatmask(index(u))=1;
    end
    tmpmask = reshape(flatmask,[Nbins,Nbins,Nslices]);
    %tmpmask = reshape(flatmask,[Nbins,Nbins,49]);
    mask = tmpmask;
    %counter = 1;
    %for u = (step-1):(step+47)
    %  mask(:,:,u) = tmpmask(:,:,counter);
    %  counter = counter+1;
    %end
    % Write segmentation result (-> mask) to file
    name = strcat('mask_', filename,'_',num2str(i),'.raw');
    fid = fopen(name, 'w');
    fwrite(fid,mask,'float32');
    fclose(fid);
  else  
    % load sino-mask and directly apply to sino
%    Mask = zeros(Nbins,Nproj,Nslices,'double');
%    name = strcat('mask_proj_01Blank',num2str(i),'.raw');
%    fid  = fopen(name,'r');
%    for l=1:Nslices
%      Mask2D = fread(fid,[Nbins,Nproj],'float32');
%      Mask(:,:,l)  = Mask(:,:,l)  + Mask2D;
%    end
%    fclose(fid);
%    SSRBSino = SSRBSino.*Mask;
%    name = strcat('SSRB_seg_',filename,'_',num2str(i),'.raw');
%    fid = fopen(name,'w');
%    fwrite(fid,SSRBSino,'float32');
%    fclose(fid);
  end
end

if strcmp(filename, '01Blank')
  % Reproject the Segmentation to get Sinogram with reduced scatter/randoms
  PETprojection(filename,Nparts);
end

% Add up all parts to get Sinogram of whole acquisition
SSRBSino = zeros(Nbins,Nproj,Nslices,'double');
%MaskSum  = zeros(Nbins,Nbins,Nslices,'double');
for i=1:Nparts
  name1 = strcat('SSRB_cor_',filename,'_',num2str(i),'.raw');
  %name2 = strcat('mask_',filename,'_',num2str(i),'.raw');
  fid1  = fopen(name1,'r');
  %fid2  = fopen(name2,'r');
  for l=1:Nslices
    Sino2D = fread(fid1,[Nbins,Nproj],'float32');
    SIN2D  = double(Sino2D);
    %Mask2D = fread(fid2,[Nbins,Nbins],'float32');
    for j=1:Nbins
      for k=1:Nproj
        if isnan(SIN2D(j,k));
          SIN2D(j,k) = 0.0;
        end
      end
    end
    SSRBSino(:,:,l) = SSRBSino(:,:,l) + SIN2D;
    %MaskSum(:,:,l)  = MaskSum(:,:,l)  + Mask2D;
  end
  fclose(fid1);
  %fclose(fid2);
end
% Write Sinogram to file
name = strcat('SSRB_complete_', filename,'.raw');
fid3 = fopen(name, 'w');
fwrite(fid3,SSRBSino,'float32');
fclose(fid3);
%name = strcat('Mask_complete_', filename,'.raw');
%fid4 = fopen(name, 'w');
%fwrite(fid4,MaskSum,'float32');
%fclose(fid4);

end

function saveMask(mask,filename)
mask = logical(mask);

fid=fopen(filename,'w');
fwrite(fid,mask,'ubit1');
fclose(fid);

% To Read it again use:
%fid=fopen(filename,'r');
%mask=fread(fid,[a b],'ubit1=>double');
%mask=fread(fid,[a b],'ubit1=>logical');
end
