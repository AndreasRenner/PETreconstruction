function SSRB_mMR_add(filename, Nscans)

% Values in [s] of time from start of total scan to middle of
% respective scan
corBlankFast = [98,218.5,340,464.5,583,699,808.5,926,1047,1156.5];
corTransPhantom1 = [60,185];
corTransPhantom2 = [80,197];
corTransRQ = [71.5,185.5];

% Basic Parameters of Siemens Biograph mMR
Nbins   = 344;         % Number of radial bins (NRAD)
Nproj   = 252;         % Number of projections (NANG)
Nslices = 127;         % Number of slices (2*Nrings -1)

% SSRB Algorithm: creates a rebinned direct sinogram  
SSRBSino = zeros(Nbins,Nproj,Nslices,'double');

% Add time from filling of Pellet to start of total scan
if     strcmp(filename, '05BlankFast')
  decayCor = corBlankFast + 3471;
elseif strcmp(filename, '02TransPhantom1')
  decayCor = corTransPhantom1 + 2234;
elseif strcmp(filename, '03TransPhantom2')
  decayCor = corTransPhantom2 + 2630;
elseif strcmp(filename, '04TransRQ')
  decayCor = corTransRQ + 3039;
else
  decayCor = zeros(Nscans,1);
end

halflifeFDG = 6586.2;   % halflife of FDG in [s]
decayF = zeros(length(decayCor),1);

for i=1:Nscans
  decayF(i) = 2^(-decayCor(i)/halflifeFDG);
  fprintf('Decay Faktor for Scan %u is: %u\r',i,decayF(i));
  name = strcat('sino_SSRB_',filename,num2str(i),'.raw');
  fid  = fopen(name,'r');
  for j=1:Nslices
    Sino2D = fread(fid,[Nbins,Nproj],'float32');
    SIN2D  = double(Sino2D)/decayF(i);
    SSRBSino(:,:,j) = SSRBSino(:,:,j) + SIN2D;
  end
  fclose(fid);
end

name = strcat('sino_SSRB_decayCor_', filename, '.raw');
fid = fopen(name,'w');
fwrite(fid,SSRBSino,'float32');
fclose(fid);

end

