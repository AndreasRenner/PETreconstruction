function CTuMap=carneyBilinearScaling(pathOfCTdicom)

% Parameters for Bilinear scaling, choose a, b and BP.
cd(pathOfCTdicom);
CTfiles=dir;
CTfiles=CTfiles(arrayfun(@(x) x.name(1), CTfiles) ~= '.'); % reading the files inside the folder.

% Load the dicom images in a volume
[volume,~,~]=dicomreadVolume(cd);
CTimg(:,:,:)=volume(:,:,1,:);

% Get the information about important tags from dicom-header
CTdcmInfo=dicominfo(CTfiles(1).name);
TubeVoltage=CTdcmInfo.KVP;
PixelSpacing=CTdcmInfo.PixelSpacing;
PixelSpacing(3)=CTdcmInfo.SliceThickness;
%StudyTime=CTdcmInfo.StudyTime
%FileModDate=CTdcmInfo.FileModDate
if isfield(CTdcmInfo,'RescaleSlope')
  RescaleSlope=CTdcmInfo.RescaleSlope;
else
  fprintf('Field does not exist!\r');
  RescaleSlope=1;
end
if isfield(CTdcmInfo,'RescaleIntercept')
  RescaleIntercept=CTdcmInfo.RescaleIntercept;
else
  fprintf('Field does not exist!\r');
  RescaleIntercept=0;
end
ImageOrientation=CTdcmInfo.ImageOrientationPatient;
kvector = zeros(3,1);
for i=1:3
  ivector(i) = ImageOrientation(i);
  jvector(i) = ImageOrientation(3+i);
  if ivector(i) || jvector(i)
    kvector(i) = 0;
  else
    kvector(i) = 1;
  end
end
transA = zeros(4,4);
transA(4,4) = 1;
for i=1:3
  transA(i,2) = ivector(i);
  transA(i,1) = jvector(i);
  transA(i,3) = kvector(i);
end

transA = affine3d(transA);


% Convert Dicom Volume in 3D Volume with correct units
if RescaleSlope==1 && RescaleIntercept==0
  CTimg(:,:,:)=volume(:,:,1,:);
else
  CTimg(:,:,:)=RescaleSlope*double(volume(:,:,1,:))+RescaleIntercept;
end


%CTimgToScale = CTimg;
% Transform Image to Patient Orientation
CTimgToScale = imwarp(CTimg,transA);

switch TubeVoltage % these values are obtained from the literature (Carney et.al, 2006, Transforming CT images for attenuation correction in PET/CT
  case 80
    a  = 3.64 * 10^-5;
    b  = 6.26 * 10^-2;
    BP = 1050; % (HU+1000) based on scaled components (0 to 4000)
    disp(a);
    disp(b);
    disp(BP);
  case 100
    a  = 4.43 * 10^-5;
    b  = 5.44 * 10^-2;
    BP = 1052; % (HU+1000) based on scaled components (0 to 4000)
    disp(a);
    disp(b);
    disp(BP);
  case 110
    a  = 4.92 * 10^-5;
    b  = 4.88 * 10^-2;
    BP = 1043; % (HU+1000) based on scaled components (0 to 4000)
    disp(a);
    disp(b);
    disp(BP);
  case 120
    a  = 5.10 * 10^-5;
    b  = 4.71 * 10^-2;
    BP = 1070; % (HU+1000) based on scaled components (0 to 4000)
    % old number was 1047; according to Clays its 1070
    disp(a);
    disp(b);
    disp(BP);
  case 130
    a  = 5.51 * 10^-5;
    b  = 4.24 * 10^-2;
    BP = 1037; % (HU+1000) based on scaled components (0 to 4000)
    disp(a);
    disp(b);
    disp(BP);
  case 140
    a  = 5.64 * 10^-5;
    b  = 4.08 * 10^-2;
    BP = 1030; % (HU+1000) based on scaled components (0 to 4000)
    disp(a);
    disp(b);
    disp(BP);
end

% Scaling the CT using the bilinear coefficients.
CTbelowBP=double(CTimgToScale).*double(CTimgToScale<=(BP));
CTbelowBPto511KeV=(9.6*10^-5).*CTbelowBP.*double(CTimgToScale<=(BP)); % values in LAC, units cm^-1
CTaboveBP=double(CTimgToScale).*double(CTimgToScale>(BP));
CTaboveBPto511KeV=(a.*CTaboveBP+b).*double(CTimgToScale>(BP)); % values in LAC, units cm^-1
CTto511=CTbelowBPto511KeV+CTaboveBPto511KeV;
%CTto511=round(10000.*CTto511); % specially scaled for siemens mMR reconstruction. 
CTto511(CTto511<0)=0;

%% segment bed
% start by dividing the image in foreground and background
levelThresh = multithresh(CTto511);
segmentation = imquantize(CTto511,levelThresh); %0.088 for Pig
segmentation = segmentation-1;
for i=1:size(segmentation,3)
  segmentation(:,:,i)=imfill  (segmentation(:,:,i),'holes');
  segmentation(:,:,i)=imerode (segmentation(:,:,i),strel('disk',5));
  segmentation(:,:,i)=imdilate(segmentation(:,:,i),strel('disk',6));
end
% search for the largest volume
conncomp=bwconncomp(segmentation);
[~,maxcell]=max(cellfun(@numel,conncomp.PixelIdxList));
% create a mask with largest volume
mask=zeros(size(segmentation));
mask(conncomp.PixelIdxList{1,maxcell})=1;
% apply mask to CTuMap
CTuMap=CTto511.*mask;
% 
CTuMap = CTuMap*10000;

% To convert MRAC to AC Map
%CTuMap = double(CTimgToScale);%/10000;

%% write result to file
name = 'Pig2_UTEnew.raw';
fid = fopen(name,'w');
%fwrite(fid,CTuMap,'float32');
fwrite(fid,UTEnew,'float32');
%fwrite(fid,uMap,'float32');
fclose(fid);


%%
% Other stuff; playing around; for evaluation of phantom study

%% read old results
%name = 'ConvertedCTtoACMapSmooth.raw';
name = 'Pig2_ECAT.raw';
Nbins=126;
Nproj=192;
Nslices=128;
uMap = zeros(Nbins,Nproj,Nslices,'double');
fid = fopen(name,'r');
for l=1:Nslices
  uMap2D = fread(fid,[Nbins,Nproj],'float32');
  UMAP2D = double(uMap2D);
  uMap(:,:,l) = uMap(:,:,l) + UMAP2D;
end
fclose(fid);

dicomImages=dir;
dicomImages=dicomImages(~[dicomImages.isdir]); %removes the variables which do not have a name.

for lp=1:length(dicomImages)
    info = dicominfo(dicomImages(lp).name); % read the whole file size of a single size
    rowsInfo=double(info.Rows);columnsInfo=double(info.Columns);
    ImpInfo=ReadDicomElementList(dicomImages(lp).name); % store the structs which are generated by dicomdisp
    fid=fopen(dicomImages(lp).name,'r+'); % open the file for reading
    OffSetData=info.FileSize-ImpInfo(end).location; 
    status=fseek(fid,-(OffSetData+(rowsInfo*columnsInfo*2)),'eof'); % always go from 'end of the file', place the cursor where the pixel data starts, including the offset
    % find out where the pixel data really ends, works for siemens.
    A=uMap(:,:,info.InstanceNumber);
    A = A.';
    count=fwrite(fid,A(:),'uint16'); % write the CT uMap to the MR uMap slices.
    fclose(fid); % always close the file.
end

%% push data to DICOM
% replaces the image-data of the original MRAC Dicom with the wanted
% image-data from an other modality
umap = UTEumap;
dicomImages=dir;
dicomImages=dicomImages(~[dicomImages.isdir]); %removes the variables which do not have a name.
%umap = rot90(umap,2);
%umap = flipdim(umap,1);
for lp=1:length(dicomImages)
    info = dicominfo(dicomImages(lp).name); % read the whole file size of a single size
    rowsInfo = double(info.Rows);
    columnsInfo = double(info.Columns);
    %ImpInfo = ReadDicomElementList(dicomImages(lp).name); % store the structs which are generated by dicomdisp
    fid=fopen(dicomImages(lp).name,'r+'); % open the file for reading
    %OffSetData=info.FileSize-ImpInfo(end).location; 
    status=fseek(fid,-(rowsInfo*columnsInfo*2),'eof'); % always go from 'end of the file', place the cursor where the pixel data starts, including the offset
    % find out where the pixel data really ends, works for siemens.
    A = umap(:,:,info.InstanceNumber); % get the corresponding CT uMap slice based on the instance number.
    A = A.';
    count=fwrite(fid,A(:),'uint16'); % write the CT uMap to the MR uMap slices.
    fclose(fid); % always close the file.
end


end
