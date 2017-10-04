function main(filenameB,filenameT)

% Basic Parameters of Siemens Biograph mMR
Nbins   = 344;         % Number of radial bins (NRAD)
Nproj   = 252;         % Number of projections (NANG)
Nslices =  13;         % Number of SSRB sinogram planes
Nparts  = 500;         % Number of parts per scan

% Load information about Blank acquisition
if     strcmp(filenameB,'05BlankFast')
  decayCor = 3471;
  NscansB  = 10;
  faktor   = [1.05,1.05,1.05,1,0.98,0.9,1.08,0.9,1,1];
elseif strcmp(filenameB,'01Blank')
  decayCor = 930;
  NscansB  = 16;
  faktor   = [1.16,1.14,1.11,1.08,1.07,1.05,1.06,1.01,1.02,1,0.94,0.92,0.92,0.9,0.88];
  direction= 1;
end

sinoTotal = zeros(Nparts,Nbins,Nproj,Nslices);

% For first Scan do not perform adaption of position
% -> first Scan is used as reference
[sino,index] = readScan(filenameB,1,NscansB,Nparts,faktor(1),0);
sinoTotal = sinoTotal + sino;
for i=1:Nparts
  fprintf('For part %i max index is %i\r', i, index(i));
end

[line1,line2] = sino2line(sinoTotal,Nparts);

% For other Scans perform adaption of position
kumFaktor = faktor(1);
for i=2:(NscansB-1)
  offset = uint64(ceil(Npos*kumFaktor))*4;
  [sino,index] = readScan(filenameB,i,NscansB,Nparts,faktor(i),offset);
  kumFaktor = kumFaktor + faktor(i);
end

% Load information about Transmission acquisition
if     strcmp(filenameT,'02TransPhantom1')
  decayCor = 2234;
  NscansT  = 2;
elseif strcmp(filenameT,'03TransPhantom2')
  decayCor = 2630;
  NscansT  = 2;
elseif strcmp(filenameT,'04TransRQ')
  decayCor = 3039;
  NscansT  = 2;
elseif strcmp(filenameT,'02TransPhantom2')
  decayCor = 2910;
  NscansT  = 6;
  faktor   = [1.05,1.03,0.97,0.98,0.95,0.95];
  direction= 1;
elseif strcmp(filename, '03TransPhantom1')
  decayCor = 3668;
  NscansT  = 7;
  faktor   = [1.05,1.05,1.02,1,0.98,0.96,1];
  direction=1;
elseif strcmp(filename, '04TransPhantom2HOT')
  decayCor = 6231;
  NscansT  = 7;
  faktor   = [0.9,0.9,1.5,1.42,0.82,0.75,1];
  direction= 0;
end

for i=1:NscansT
    
  if direction
    direction = 0;
  else
    direction = 1;
  end
end

end