function main(filenameB,filenameT)

%% Basic Parameters
% of Siemens Biograph mMR
Nbins   = 344;         % Number of radial bins (NRAD)
Nproj   = 252;         % Number of projections (NANG)
Nslices =  13;         % Number of SSRB sinogram planes
% other parameters
Nparts  = 500;         % Number of parts per scan
halflifeFDG = 6586.2;   % halflife of FDG in [s]


%% B L A N K   S C A N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load information about Blank acquisition
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

%% Allocate memory for Blank acquisition
decayF     = zeros(Nparts,1);
sinoTotalB = zeros(Nparts,Nbins,Nproj,Nslices);
tstartB    = zeros(Nparts,NscansB);
tstopB     = zeros(Nparts,NscansB);

%% Read First Scan of Blank acquisition
% For first Scan do not perform adaption of position
% -> first Scan is used as reference
[sino,indexref,tstart,tstop] = readScan(filenameB,1,NscansB,Nparts, ...
    faktor(1),0,direction,0,0);

%% DECAY CORRECTION
durrationRef = tstop-tstart;
for i=1:Nparts
  decayF(i) = 2^(-(decayCor+tstart(i)+durrationRef(i)/2)/halflifeFDG);
  sino(i,:,:,:) = sino(i,:,:,:)*decayF(i);
end
% No scan-time correction for this scan, as this scan-time
% is used as reference and a constant value for all parts

%% Save files for visual check
%for i=1:Nparts
%  name = strcat('Scan1_',num2str(i));
%  SSRBSino(:,:,:)=sino(i,:,:,:);
%  newName = strcat('sino_SSRB_', name, '.raw');
%  fid = fopen(newName,'w');
%  fwrite(fid,SSRBSino,'float32');
%  fclose(fid);
%  clear SSRBSino;
%end

%%
sinoTotalB = sinoTotalB + sino;

%%
tstartB(:,1) = tstart(:);
tstopB(:,1)  = tstop(:);

%% Calculate sinogram window times for reference
tic
[line1and2] = sino2line(sino,Nparts,1);
toc

%% For other Scans perform adaption of position
kumFaktor = faktor(1);
for i=2:(NscansB-1)
  if direction
    direction = 0;
  else
    direction = 1;
  end
  [sino,~,tstart,tstop] = readScan(filenameB,i,NscansB,Nparts, ...
      faktor(i),kumFaktor,direction,line1and2,indexref);
  
  %% DECAY CORRECTION + SCAN-TIME CORRECTION
  durration = tstop-tstart;
  scanTcor  = durrationRef./durration;
  for j=1:Nparts
    decayF(j) = 2^(-(decayCor+tstart(j)+durration(j)/2)/halflifeFDG);
    sino(j,:,:,:) = sino(j,:,:,:)*decayF(j)*scanTcor(j);
  end
  
  %% Save files for visual check
  for j=1:Nparts
    name = strcat('Scan',num2str(i),'_',num2str(j));
    SSRBSino(:,:,:)=sino(j,:,:,:);
    newName = strcat('sino_SSRB_', name, '.raw');
    fid = fopen(newName,'w');
    fwrite(fid,SSRBSino,'float32');
    fclose(fid);
    clear SSRBSino;
  end
  
  %% add new sinograms to total sinogram
  sinoTotalB = sinoTotalB + sino;
  %% add times for decay correction
  tstartB(:,i) = tstart(:);
  tstopB(:,i)  = tstop(:);
  %%
  kumFaktor = kumFaktor + faktor(i);
end


%% T R A N S M I S S I O N   S C A N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load information about Transmission acquisition
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

%% Allocate memory for Transmission acquisition
decayF     = zeros(Nparts,1);
sinoTotalT = zeros(Nparts,Nbins,Nproj,Nslices);
tstartT    = zeros(Nparts,NscansB);
tstopT     = zeros(Nparts,NscansB);

%% For transmission scans perform adaption of position
kumFaktor = 0;
for i=1:NscansT
    
  [sino,~,tstart,tstop] = readScan(filenameT,i,NscansT,Nparts, ...
      faktor(i),kumFaktor,direction,line1and2,indexref);
  
  if direction
    direction = 0;
  else
    direction = 1;
  end
  
  %% DECAY CORRECTION + SCAN-TIME CORRECTION
  durration = tstop-tstart;
  scanTcor  = durrationRef./durration;
  for j=1:Nparts
    decayF(j) = 2^(-(decayCor+tstart(j)+durration(j)/2)/halflifeFDG);
    sino(j,:,:,:) = sino(j,:,:,:)*decayF(j)*scanTcor(j);
  end
  
  %% Save files for visual check
  for j=1:Nparts
    name = strcat('ScanTrans',num2str(i),'_',num2str(j));
    SSRBSino(:,:,:)=sino(j,:,:,:);
    newName = strcat('sino_SSRB_', name, '.raw');
    fid = fopen(newName,'w');
    fwrite(fid,SSRBSino,'float32');
    fclose(fid);
    clear SSRBSino;
  end
  
  %% add new sinograms to total sinogram
  sinoTotalT = sinoTotalT + sino;
  %% add times for decay correction
  tstartT(:,i) = tstart(:);
  tstopT(:,i)  = tstop(:);

  %%
  kumFaktor = kumFaktor + faktor(i);
end

end