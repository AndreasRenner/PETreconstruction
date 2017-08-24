function SSRB_mMR_add(filename)

% D E C A Y   C O R R E C T I O N
% Values in [s] of time from start of total scan to middle of
% respective scan
corBlankFast = [98,218.5,340,464.5,583,699,808.5,926,1047,1156.5];
corTransPhantom1 = [60,185];
corTransPhantom2 = [80,197];
corTransRQ = [71.5,185.5];
cor01Blank = [67.8,173.7,372.9,481.5,589.8,693.5,801.5,910,1024.3,1134.8,1246.5,1363.4,1466.3,1572.2,1674];
cor02Trans = [67.1,171.4,275.2,376.4,479.4,624.9];
cor03Trans = [59.9,162.8,264.3,366.4,468.3,573.9,676.7];
% S C A N T I M E   C O R R E C T I O N
% Values in [ms] (calculated with tstop-tstart)
scanTime01Blank = [80170,80310,79230,80030,79860,81460,80720,80980,79060,79710,79130,79680,79080,79520,77570];
scanTime02Trans = [79310,80270,78470,79510,78580,79040];
scanTime03Trans = [78710,79360,78430,79370,78750,79440,79740];

% Basic Parameters of Siemens Biograph mMR
Nbins   = 344;         % Number of radial bins (NRAD)
Nproj   = 252;         % Number of projections (NANG)
Nslices = 127;         % Number of slices (2*Nrings -1)

% SSRB Algorithm: creates a rebinned direct sinogram  
SSRBSino = zeros(Nbins,Nproj,Nslices,'double');

% Add time from filling of Pellet to start of total scan
if     strcmp(filename, '05BlankFast')
  decayCor = corBlankFast + 3471;
  Nscans = 10;
elseif strcmp(filename, '02TransPhantom1')
  decayCor = corTransPhantom1 + 2234;
  Nscans = 2;
elseif strcmp(filename, '03TransPhantom2')
  decayCor = corTransPhantom2 + 2630;
  Nscans = 2;
elseif strcmp(filename, '04TransRQ')
  decayCor = corTransRQ + 3039;
  Nscans = 2;
elseif strcmp(filename, '01Blank')
  decayCor = cor01Blank + 930;
  Nscans = 15;
  scanTime = scanTime01Blank;
elseif strcmp(filename, '02TransPhantom2')
  decayCor = cor02Trans + 2910;
  Nscans = 6;
  scanTime = scanTime02Trans;
elseif strcmp(filename, '03TransPhantom1')
  decayCor = cor03Trans + 3668;
  Nscans = 7;
  scanTime = scanTime03Trans;
else
  Nscans = 10;
  decayCor = zeros(Nscans,1);
end

halflifeFDG = 6586.2;   % halflife of FDG in [s]

for i=1:Nscans
  decayF    = 1/2^(-decayCor(i)/halflifeFDG);
  scanTimeF = 100000/scanTime(i);
  fprintf('Decay Faktor for Scan %u is: %u\r',i,decayF);
  name=strcat('sino_SSRB_',filename,num2str(i),'randomSubstracted_noGapFilling.raw');
  fid =fopen(name,'r');
  for j=1:Nslices
    Sino2D = fread(fid,[Nbins,Nproj],'float32');
    % Apply Scantime-Cor + Decay-Cor + Scannumber-Cor
    SIN2D  = double(Sino2D)*scanTimeF*decayF/Nscans;
    SSRBSino(:,:,j) = SSRBSino(:,:,j) + SIN2D;
  end
  fclose(fid);
end

% Read Mask for Sinogram and Mask for Ratio
SinoMaskS= zeros(Nbins,Nproj,Nslices,'double');
nameMaskS = strcat('mask_Tube_2017_08.raw');
fidMaskS  = fopen(nameMaskS,'r');
for i=1:Nslices
  SinoMaskS(:,:,i) = fread(fidMaskS,[Nbins,Nproj],'float32');
end
fclose(fidMaskS);

% Apply Mask
SSRBSino = SSRBSino.*SinoMaskS;

name = strcat('sino_SSRB_seg_Cor_noGap_',num2str(Nscans),filename,'.raw');
fid = fopen(name,'w');
fwrite(fid,SSRBSino,'float32');
fclose(fid);

end

