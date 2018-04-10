function main(filenameB,filenameT)

%% Basic Parameters
% of Siemens Biograph mMR
Nbins   = 344;         % Number of radial bins (NRAD)
Nproj   = 252;         % Number of projections (NANG)
Nslices =  13;         % Number of SSRB sinogram planes around pellet
NSliTot = 127;         % Total Number of SSRB sinogram planes
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
  NscansB  = 16; % last Scan is not comlete
  faktor   = [1.16,1.14,1.11,1.08,1.07,1.05,1.06,1.01,1.02,1,0.94,0.92,0.92,0.9,0.88];
  direction= 1;
elseif strcmp(filenameB,'blankPig1')
  decayCor = 0;
  NscansB  = 22; % last Scan is not complete
  faktor   = [1.14,1.14,1.12,1.09,1.08,1.08,1.06,1.08,1.03,1.03, ...
      1.03,1,1,0.98,0.96,0.96,0.94,0.96,0.92,0.93,0.9,1,1];
  direction= 0;
elseif strcmp(filenameB,'blankPig2')
  decayCor = 0;
  NscansB  = 22;
  faktor   = [1.11,1.12,1.1,1.07,1.08,1.04,1.03,1.03,1.04,1,1, ...
      0.99,0.99,0.98,0.97,0.94,0.94,0.93,0.92,0.93,0.9,1,1];
  direction= 0;
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
    faktor(1),0,direction,0,0,0);

%% DECAY CORRECTION
durrationRef = tstop-tstart;
for i=1:Nparts
  decayF(i) = 2^(-(decayCor+(tstart(i)+durrationRef(i)/2)/1000)/halflifeFDG);
  sino(i,:,:,:) = sino(i,:,:,:)/decayF(i);
end
% No scan-time correction for this scan, as this scan-time
% is used as reference and a constant value for all parts

%% Save files for visual check
cd('/media/andreas/CorneaOCT_StudyData2/Backup/Leipzig/Evaluation/ScansBlankPig1/sinoScan1');
for i=1:Nparts
  name = strcat('Scan1_',num2str(i));
  SSRBSino(:,:,:)=sino(i,:,:,:);
  newName = strcat('sino_SSRB_', name, '.raw');
  fid = fopen(newName,'w');
  fwrite(fid,SSRBSino,'float32');
  fclose(fid);
  clear SSRBSino;
end
cd('/home/andreas/code/PETreconstruction');

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
      faktor(i),kumFaktor,direction,dirRef,line1and2,indexref);
  
  %% DECAY CORRECTION + SCAN-TIME CORRECTION
  durration = tstop-tstart;
  scanTcor  = durrationRef./durration;
  for j=1:Nparts
    decayF(j) = 2^(-(decayCor+(tstart(j)+durration(j)/2)/1000)/halflifeFDG);
    sino(j,:,:,:) = sino(j,:,:,:)/decayF(j)*scanTcor(j);
  end
  
  %% Save files for visual check
  pathname = strcat('/media/andreas/CorneaOCT_StudyData2/Backup/Leipzig/Evaluation/ScansBlankPig1/sinoScan',num2str(i));
  try
    cd(pathname);
  catch ME
    ME;
  end
  for j=1:Nparts
    name = strcat('Scan',num2str(i),'_',num2str(j));
    SSRBSino(:,:,:)=sino(j,:,:,:);
    newName = strcat('sino_SSRB_', name, '.raw');
    fid = fopen(newName,'w');
    fwrite(fid,SSRBSino,'float32');
    fclose(fid);
    clear SSRBSino;
  end
  cd('/home/andreas/code/PETreconstruction');
  
  %% add new sinograms to total sinogram
  sinoTotalB = sinoTotalB + sino;
  %% add times for decay correction
  tstartB(:,i) = tstart(:);
  tstopB(:,i)  = tstop(:);
  %%
  kumFaktor = kumFaktor + faktor(i);
end

%% Read total Sinogram of Blank Scan from single Scan files
sinoTotalB = zeros(Nparts,Nbins,Nproj,Nslices);
cd('/media/andreas/CorneaOCT_StudyData2/Backup/Leipzig/Evaluation/ScansBlankPig2');
for i=1:NscansB
  name = strcat('sinoScan',num2str(i));
  cd(name)
  sinoTotalTmp = zeros(Nparts,Nbins,Nproj,Nslices);
  for j=1:Nparts
    sinoname = strcat('sino_SSRB_Scan',num2str(i),'_',num2str(j),'.raw');
    SSRBSino = zeros(Nbins,Nproj,Nslices,'double');
    fid = fopen(sinoname,'r');
    for l=1:Nslices
      Sino2D = fread(fid,[Nbins,Nproj],'float32');
      SIN2D = double(Sino2D);
      SSRBSino(:,:,l) = SSRBSino(:,:,l) + SIN2D;
    end
    fclose(fid);
    sinoTotalTmp(j,:,:,:) = SSRBSino(:,:,:);
  end
  sinoTotalB = sinoTotalB + sinoTotalTmp;
  cd ..
end
clear sinoTotalTmp;
cd('/home/andreas/code/PETreconstruction');

%% Correct for number of scans
sinoTotalB = sinoTotalB./i;

%% Write total Sinogram of Blank Scan to file
for j=1:Nparts
  name = strcat(filenameB,'_',num2str(j));
  SSRBSino(:,:,:)=sinoTotalB(j,:,:,:);
  newName = strcat('sino_SSRB_', name, '.raw');
  fid = fopen(newName,'w');
  fwrite(fid,SSRBSino,'float32');
  fclose(fid);
  clear SSRBSino;
end

%% Read total Sinogram of Blank Scan from file
sinoTotalB = zeros(Nparts,Nbins,Nproj,Nslices);
for j=1:Nparts
  SSRBSino = zeros(Nbins,Nproj,Nslices,'double');
  name = strcat(filenameB,'_',num2str(j));
  filename = strcat('sino_SSRB_', name, '.raw');
  fid = fopen(filename,'r');
  for i=1:Nslices
    Sino2D = fread(fid,[Nbins,Nproj],'float32');
    SIN2D = double(Sino2D);
    SSRBSino(:,:,i) = SSRBSino(:,:,i) + SIN2D;
  end
  fclose(fid);
  sinoTotalB(j,:,:,:) = SSRBSino(:,:,:);
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
elseif strcmp(filenameT, '03TransPhantom1')
  decayCor = 3668;
  NscansT  = 7;
  faktor   = [1.05,1.05,1.02,1,0.98,0.96,1];
  direction= 1;
elseif strcmp(filenameT, '04TransPhantom2HOT')
  decayCor = 6231;
  NscansT  = 7;
  faktor   = [0.9,0.9,1.5,1.42,0.82,0.75,1];
  direction= 0;
elseif strcmp(filenameT, 'transPig1')
  decayCor = 3193;
  NscansT  = 12;
  faktor   = [1.2,1.1,1.1,1,1,1,1,0.98,0.92,0.88,0.9,1,1];
  direction= 0;
elseif strcmp(filenameT, 'transPig2')
  decayCor = 3672;
  NscansT  = 12;
  faktor   = [1.1,1.11,1.04,1.06,1.02,1,0.98,0.96,0.94,0.92,0.9,1,1];
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
      faktor(i),kumFaktor,direction,dirRef,line1and2,indexref);
  
  if direction
    direction = 0;
  else
    direction = 1;
  end
  
  %% DECAY CORRECTION + SCAN-TIME CORRECTION
  durration = tstop-tstart;
  scanTcor  = durrationRef./durration;
  for j=1:Nparts
    decayF(j) = 2^(-(decayCor+(tstart(j)+durration(j)/2)/1000)/halflifeFDG);
    sino(j,:,:,:) = sino(j,:,:,:)/decayF(j)*scanTcor(j);
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

%% Correct for number of scans
sinoTotalT = sinoTotalT./i;

%% Write total Sinogram of Transmission Scan to file
for j=1:Nparts
  name = strcat(filenameT,'_',num2str(j));
  SSRBSino(:,:,:)=sinoTotalT(j,:,:,:);
  newName = strcat('sino_SSRB_', name, '.raw');
  fid = fopen(newName,'w');
  fwrite(fid,SSRBSino,'float32');
  fclose(fid);
  clear SSRBSino;
end

%% Read total Sinogram of Trans Scan from file
sinoTotalT = zeros(Nparts,Nbins,Nproj,Nslices);
for j=1:Nparts
  SSRBSino = zeros(Nbins,Nproj,Nslices,'double');
  name = strcat(filenameT,'_',num2str(j));
  filename = strcat('sino_SSRB_', name, '.raw');
  fid = fopen(filename,'r');
  for i=1:Nslices
    Sino2D = fread(fid,[Nbins,Nproj],'float32');
    SIN2D = double(Sino2D);
    SSRBSino(:,:,i) = SSRBSino(:,:,i) + SIN2D;
  end
  fclose(fid);
  sinoTotalT(j,:,:,:) = SSRBSino(:,:,:);
end

%% Read Emission Data - Multiple Emission Sinograms
NpartE  = 13;
%tstart  = [ 0,141,283,423,563,703,844, 982,1123,1263,1406,1546,1688]; %TransPig2
%tstop   = [64,205,347,487,627,767,908,1046,1187,1328,1470,1610,1799];
tstart  = [0,146,299,444,591,733,876,1021,1162,1302,1443,1586,1737]; %TransPig1
tstop   =[63,209,362,507,654,796,940,1084,1225,1365,1507,1649,1799];
acqtime = tstop-tstart; %for each sinogram
%
sinoEm  = zeros(Nbins,Nproj,NSliTot);
for i=1:NpartE
  % (1) DECAY CORRECTION
  decayF(i) = 2^(-(decayCor+tstart(i)+acqtime(i)/2)/halflifeFDG);
  % (2) SCAN-TIME CORRECTION
  % meantime could be replaced by the actual time of each
  % transmission scan
  meantime = 0.131; % BlankPig1
  %meantime = 0.131; % BlankPig2
  timeF(i) = meantime/acqtime(i);
  % (3) SCAN-NUMBER CORRECTION
  totalCorF(i) = timeF(i)/decayF(i)/NpartE;
  nameEm = strcat('sino_SSRB_',filenameT,num2str(i),'emission.raw');
  fid = fopen(nameEm,'r');
  for j=1:NSliTot
    sinoEm2D = fread(fid,[Nbins,Nproj],'float32');
    SIN2D = double(sinoEm2D)*totalCorF(i);
    sinoEm(:,:,j) = sinoEm(:,:,j) + SIN2D;
  end
  fclose(fid);
end

%% Read Emission Data - Single Emission Sinogram
tstart = 311.96;
acqtime= 317.977;
% (1) DECAY CORRECTION
decayF   = 2^(-(decayCor+tstart+acqtime/2)/halflifeFDG);
% (2) SCAN-TIME CORRECTION
% meantime could be replaced by the actual time of each
% transmission scan
%meantime = 0.16; % 01Blank
meantime = 0.131; % BlankPig2
timeF    = meantime/acqtime;
% (3) SCAN-NUMBER CORRECTION is not necessary (counts are already
%     divided by number of scans)
totalCorF = timeF/decayF;
sinoEm = zeros(Nbins,Nproj,NSliTot);
nameEm = 'sino_SSRB_04TransPhantom2HOTstart_311960_acqtime_317977.raw';
fid = fopen(nameEm,'r');
for i=1:NSliTot
  sinoEm2D = fread(fid,[Nbins,Nproj],'float32');
  SIN2D = double(sinoEm2D)*totalCorF;
  sinoEm(:,:,i) = sinoEm(:,:,i) + SIN2D;
end
fclose(fid);

%% make sino Total Emission
% this is needed for scatter reduction using line1and2
sinoTotalE = zeros(Nparts,Nbins,Nproj,Nslices);
currentEmiss = zeros(Nbins,Nproj,Nslices);
for j=1:Nparts
  l = 1;
  for i=(indexref(j)-6):(indexref(j)+6)
    if i>0 && i<=NSliTot
      % substract emission sinogram from current part
      currentEmiss(:,:,l) = sinoEm(:,:,i);
    end
    l = l+1;
  end
  % save modified current part
  sinoTotalE(j,:,:,:) = currentEmiss(:,:,:);
end

tic
%% Reduce Scatter using line1and2 for segmentation
%problemcount = 0;
% changelow is going up to the left side
% changetop is going up to the right side
% -> first entry of the vector is k-value of startingpoint
% Define object boundaries:
% objectstart < object < objectstop
% object<= objectstop
% object>= objectstart
testMask = ones(Nparts,Nbins,Nproj,Nslices);
objectstart = 115;
objectstop  = 214;
% get data for part 1
border1old(:) = line1and2(1,1,:);
border2old(:) = line1and2(1,2,:);
changelowold  = 0;
changetopold  = 0;
for j=1:Nparts
  changelownew = 0;
  changetopnew = 0;
  % get data for part j+1
%  fprintf('Part %i: ',j);
  if j<500
    border1new(:) = line1and2(j+1,1,:);
    border2new(:) = line1and2(j+1,2,:);
  end
  lowerlimitold = zeros(Nproj,1);
  upperlimitold = zeros(Nproj,1);
  lowerlimitnew = zeros(Nproj,1);
  upperlimitnew = zeros(Nproj,1);
  for k=1:Nproj
    % reduce border to the pixel corresponding to Nproj
    % -> lowerlimit corresponds to line1
    % -> upperlimit corresponds to line2
    lowerlimitold(k) = round(border1old(k*10+1));
    upperlimitold(k) = round(border2old(k*10+1));
    lowerlimitnew(k) = round(border1new(k*10+1));
    upperlimitnew(k) = round(border2new(k*10+1));
  end
  if lowerlimitold(border1old(end-1))>165
    %
    % --------------------------------------------------------
    %% CASE 1 --- for this case line1 is top line in lower part
    %fprintf('Maximum is on right side (case 1).\r');
    % --------------------------------------------------------
    %
    %% First change lowerlimit:
    %fprintf('Changing lowerlimit...\r');
    stop = 0;
    k = 1;
    if lowerlimitold(k)<objectstart
      while lowerlimitold(k)<objectstart
        lowerlimit(k) = lowerlimitold(k);
        k = k+1;
      end
    end
    if lowerlimitold(k)<=objectstop
      i = 1;
      changetopnew(1) = k;
      i = i+1;
      while lowerlimitold(k)<=objectstop
        lowerlimit(k) = (lowerlimitold(k)+upperlimitnew(k))/2;
        changetopnew(i) = lowerlimit(k);
        i = i+1;
        k = k+1;
      end
      lowerlimit(k) = (lowerlimit(k-1)+lowerlimitold(k+1))/2;
      k = k+1;
    end
    while lowerlimitold(k)>objectstop
      lowerlimit(k) = lowerlimitold(k);
      if k==Nproj
        %fprintf('Needed to stop while loop.\r');
        stop = 1;
        break
      end
      k = k+1;
    end
    if changelowold(1) && ~stop
      if k~=changelowold(1)
        %fprintf('There is a difference in k!\r');
        k = changelowold(1);
      end
      lowerlimit(k-1) = (lowerlimitold(k-2)+changelowold(2))/2;
      for i=2:length(changelowold)
        lowerlimit(k) = changelowold(i);
        k = k+1;
      end
    end
    if ~stop
      while k<=Nproj
        lowerlimit(k) = lowerlimitold(k);
        k = k+1;
      end
    end
    %
    %% Second change upperlimit:
    %fprintf('Changing upperlimit...\r');
    stop = 0;
    k = 1;
    if upperlimitold(1)<objectstart
      while upperlimitold(k)<objectstart
        upperlimit(k) = upperlimitold(k);
        k = k+1;
      end
    end
    if changetopold
      if k~=changetopold(1)
        %fprintf('There is a difference in k!\r');
        k = changetopold(1);
      end
      for i=2:length(changetopold);
        upperlimit(k) = changetopold(i);
        k = k+1;
      end
      upperlimit(k) = (changetopold(end)+upperlimitold(k+1))/2;
      k = k+1;
    else
      while upperlimitold(k)<=objectstop
        upperlimit(k) = upperlimitold(k);
        k = k+1;
      end
    end
    while upperlimitold(k)>objectstop
      upperlimit(k) = upperlimitold(k);
      if k==Nproj
        stop = 1;
        %fprintf('Needed to stop while loop.\r');
        break;
      end
      k = k+1;
    end
    %if lowerlimitnew(Nproj)<=objectstop
    %  stop = 0;
    %end
    if ~stop
      changelownew(1) = k;
      i = 2;
      while lowerlimitnew(k)>=objectstart
        upperlimit(k) = (upperlimitold(k)+lowerlimitnew(k))/2;
        changelownew(i) = upperlimit(k);
        if k==Nproj
          stop = 1;
          %fprintf('Needed to stop while loop (point2).\r');
          break;
        end
        i = i+1;
        k = k+1;
      end
    end
    if ~stop
      upperlimit(k) = (upperlimitold(k+1)+upperlimit(k-1))/2;
      k = k+1;
      while k<=Nproj
        upperlimit(k) = upperlimitold(k);
        k = k+1;
      end
    end
  else
    %
    % --------------------------------------------------------
    %% CASE 2 --- for this case line1 is low line in lower part
    %fprintf('Maximum is on left side (case 2).\r');
    % --------------------------------------------------------
    %
    %% First change lowerlimit:
    %fprintf('Changing lowerlimit...\r');
    stop = 0;
    k = 1;
    if lowerlimitold(k)>objectstop
      while lowerlimitold(k)>objectstop
        lowerlimit(k) = lowerlimitold(k);
        k = k+1;
      end
    end
    if lowerlimitold(k)>=objectstart
      if changelowold
        if k~=changelowold(1)
          %fprintf('There is a difference in k!\r');
          k = changelowold(1);
        end
        for i=2:length(changelowold);
          lowerlimit(k) = changelowold(i);
          if k==Nproj
            %fprintf('Needed to stop for loop.\r');
            stop = 1;
            break
          end
          k = k+1;
        end
        lowerlimit(k) = (lowerlimit(k-1)+lowerlimitold(k+1))/2;
        k = k+1;
      else
        while upperlimitnew(k)>=objectstart
          lowerlimit(k) = lowerlimitold(k);
          if k==Nproj
            %fprintf('Needed to stop while loop (point1).\r');
            stop = 1;
            break
          end
          k = k+1;
        end
      end
    end
    if ~stop
      while lowerlimitold(k)<objectstart
        lowerlimit(k) = lowerlimitold(k);
        if k==Nproj
          %fprintf('Needed to stop while loop (point2).\r');
          stop = 1;
          break
        end
        k = k+1;
      end
    end
    if ~stop
      if k<Nproj
        lowerlimit(k)=(lowerlimit(k-1)+(upperlimitnew(k+1)+lowerlimitold(k+1))/2)/2;
        k = k+1;
        i = 2;
        changetopnew(1) = k;
        while lowerlimitold(k)<=objectstop
          lowerlimit(k) = (lowerlimitold(k)+upperlimitnew(k))/2;
          changetopnew(i) = lowerlimit(k);
          i = i+1;
          if k==Nproj
            %fprintf('Needed to stop while loop (point3).\r');
            stop = 1;
            break;
          end
          k = k+1;
        end
      else
        lowerlimit(k) = lowerlimitold(k);
        stop = 1;
      end
    end
    if ~stop
      if k<Nproj
        lowerlimit(k) = (lowerlimitold(k+1)+lowerlimit(k-1))/2;
        k = k+1;
        while k<=Nproj
          lowerlimit(k) = lowerlimitold(k);
          k = k+1;
        end
      else
        lowerlimit(k) = lowerlimitold(k);
      end
    end
    %
    %% Second change upperlimit:
    %fprintf('Changing upperlimit...\r');
    stop = 0;
    k = 1;
    if upperlimitold(k)>objectstop
      while upperlimitold(k)>objectstop
        upperlimit(k) = upperlimitold(k);
        k = k+1;
      end
    end
    if upperlimitold(k)>=objectstart
      i = 2;
      changelownew(1) = k;
      while upperlimitold(k)>=objectstart
        upperlimit(k) = (upperlimitold(k)+lowerlimitnew(k))/2;        
        changelownew(i) = upperlimit(k);
        i = i+1;
        k = k+1;
      end
      upperlimit(k) = (upperlimit(k-1)+upperlimitold(k+1))/2;
      k = k+1;
    end
    while upperlimitold(k)<objectstart
      upperlimit(k) = upperlimitold(k);
      if k==Nproj
        %fprintf('Needed to stop while loop.\r');
        stop = 1;
        break;
      end
      k = k+1;
    end
    if changetopold(1) && ~stop
      if k~=changetopold(1)
        %fprintf('There is a difference in k!\r');
        k = changetopold(1);
      end
      upperlimit(k-1) = (upperlimitold(k-2)+changetopold(2))/2;
      for i=2:length(changetopold)
        upperlimit(k) = changetopold(i);
        if k==Nproj
          %fprintf('Needed to stop for loop.\r');
          break
        end
        k = k+1;
      end
      while k<=Nproj
        upperlimit(k) = upperlimitold(k);
        k = k+1;
      end      
    elseif ~stop
      while k<=Nproj
        upperlimit(k) = upperlimitold(k);
        k = k+1;
      end
    end
    %
  end
  clear changelowold;
  clear changetopold;
  changelowold = changelownew;
  changetopold = changetopnew;
  clear changelownew;
  clear changetopnew;
  border1old = border1new;
  border2old = border2new;
  %if length(upperlimit)==Nproj && length(lowerlimit)==Nproj ...
  %        && min(upperlimit)>0 && min(lowerlimit)>0
  %  fprintf('DONE!\r');
  %else
  %  if length(upperlimit)~=Nproj
  %    fprintf('Problem with upperlimit! Length: %i\r',length(upperlimit));
  %  end
  %  if length(lowerlimit)~=Nproj
  %    fprintf('Problem with lowerlimit! Length: %i\r',length(lowerlimit));
  %  end
    if ~min(upperlimit)
  %    fprintf('Zeros in upperlimit!\r');
      % find zeros and fill them with upperlimitold
      emptyupper = find(upperlimit==0);
      for i=1:length(emptyupper)
        upperlimit(emptyupper(i)) = upperlimitold(emptyupper(i));
      end
    end
    if ~min(lowerlimit)
   %   fprintf('Zeros in lowerlimit!\r');
      % find zeros and fill them with lowerlimitold
      emptylower = find(lowerlimit==0);
      for i=1:length(emptylower)
        lowerlimit(emptylower(i)) = lowerlimitold(emptylower(i));
      end
    end
   % if length(upperlimit)==Nproj && length(lowerlimit)==Nproj ...
   %         && min(upperlimit)>0 && min(lowerlimit)>0
   %   fprintf('DONE!\r');
   % else
   %   fprintf('Sill a problem!\r');
   %   problemcount = problemcount + 1;
   % end
  %end
  %fprintf('\r');
  %
  for i=1:Nslices
    for k=1:Nproj
      for l=1:(round(lowerlimit(k))-1)
        sinoTotalT(j,l,k,i) = 0.0;
        sinoTotalB(j,l,k,i) = 0.0;
        %sinoTotalE(j,l,k,i) = 0.0;
        testMask  (j,l,k,i) = 0;
      end
      for l=round(upperlimit(k)):Nbins
        sinoTotalT(j,l,k,i) = 0.0;
        sinoTotalB(j,l,k,i) = 0.0;
        %sinoTotalE(j,l,k,i) = 0.0;
        testMask  (j,l,k,i) = 0;
      end
    end
  end
  clear upperlimit;
  clear lowerlimit;
end
%fprintf('Total number of problems: %i\r',problemcount);
toc

%% make complete sinogram using indexref
BlankComplete = zeros(Nbins,Nproj,NSliTot);
TransComplete = zeros(Nbins,Nproj,NSliTot);
EmissComplete = zeros(Nbins,Nproj,NSliTot);
%MaskComplete  = zeros(Nbins,Nproj,NSliTot);

for j=1:Nparts
  l = 1;
  sinoAddBlank(:,:,:) = sinoTotalB(j,:,:,:);
  sinoAddTrans(:,:,:) = sinoTotalT(j,:,:,:);
  sinoAddEmiss(:,:,:) = sinoTotalE(j,:,:,:);
%  sinoAddMask (:,:,:) = testMask  (j,:,:,:);
  for i=(indexref(j)-6):(indexref(j)+6)
    if i>0 && i<=NSliTot
      BlankComplete(:,:,i)=BlankComplete(:,:,i)+sinoAddBlank(:,:,l);
      TransComplete(:,:,i)=TransComplete(:,:,i)+sinoAddTrans(:,:,l);
      EmissComplete(:,:,i)=EmissComplete(:,:,i)+sinoAddEmiss(:,:,l);
%      MaskComplete (:,:,i)=MaskComplete (:,:,i)+sinoAddMask (:,:,l);
    end
    l = l+1;
  end
end

%tic
%% perform median filtering of complete sinogram
%for i=1:NSliTot
%  BlankComplete(:,:,i) = medfilt2(BlankComplete(:,:,i),[3 3]);
%  TransComplete(:,:,i) = medfilt2(TransComplete(:,:,i),[3 3]);
%end
%toc

%% inpaint nans to crystall gaps
for i=1:NSliTot
  Blank2D(:,:) = BlankComplete(:,:,i);
  Trans2D(:,:) = TransComplete(:,:,i);
  %Emiss2D(:,:) = EmissComplete(:,:,i);
  Mask2D(:,:)  = MaskComplete(:,:,i);
  for j=1:Nproj
    for k=1:Nbins
      if Blank2D(k,j)<=0 && Mask2D(k,j)
        Blank2D(k,j) = nan;
        Trans2D(k,j) = nan;
        %Emiss2D(k,j) = nan;
      end
    end
  end
  Blank2D = inpaint_nans(Blank2D,2);
  Trans2D = inpaint_nans(Trans2D,2);
  %Emiss2D = inpaint_nans(Emiss2D,2);
  BlankComplete(:,:,i) = Blank2D;
  TransComplete(:,:,i) = Trans2D;
  %EmissComplete(:,:,i) = Emiss2D;
end

%% write Blank Complete to file
name = strcat('SSRB_blank_noPost_',filenameB,'.raw');
fid = fopen(name,'w');
fwrite(fid,BlankComplete,'float32');
fclose(fid);

%% write Trans Complete to file
name = strcat('SSRB_trans_noPost_',filenameT,'.raw');
fid = fopen(name,'w');
fwrite(fid,TransComplete,'float32');
fclose(fid);

%% write Mask Complete to file
%name = strcat('SSRB_Mask_',filenameT,'.raw');
%fid = fopen(name,'w');
%fwrite(fid,MaskComplete,'float32');
%fclose(fid);

%% write Emiss Complete to file
name = strcat('SSRB_emiss_noPost',filenameT,'.raw');
fid = fopen(name,'w');
fwrite(fid,EmissComplete,'float32');
%fwrite(fid,sinoEm,'float32');
fclose(fid);

%% Calculate Ratio of Complete Blank and Transmission Scan
BlankComplete = smooth3(BlankComplete,'gaussian',[3 3 3],0.42466);
TransComplete = smooth3(TransComplete,'gaussian',[3 3 3],0.42466);
% medfilt3 was introduced in R2016b
%BlankComplete = medfilt3(BlankComplete,[5 5 3]);
%TransComplete = medfilt3(TransComplete,[5 5 3]);
EmissComplete = smooth3(EmissComplete,'gaussian',[3 3 3],0.42466);
%BlankComplete = smooth3(BlankComplete,'gaussian',[5 5 5],1);
%TransComplete = smooth3(TransComplete,'gaussian',[5 5 5],1);
%EmissComplete = smooth3(EmissComplete,'gaussian',[5 5 5],1);

% Substrat emission sinogram from transmission scan
TransComplete = TransComplete - EmissComplete;

SinoRatio = zeros(Nbins,Nproj,NSliTot,'double');

for u=1:Nbins
  for v=1:Nproj
    for w=1:NSliTot
      if TransComplete(u,v,w)<=0 || BlankComplete(u,v,w)<=0
        SinoRatio(u,v,w) = 0.0;
      elseif TransComplete(u,v,w)<0.1
        SinoRatio(u,v,w) = log(BlankComplete(u,v,w));
      else
        SinoRatio(u,v,w) = log(BlankComplete(u,v,w)/(TransComplete(u,v,w)));
      end
      if isnan(SinoRatio(u,v,w))
        SinoRatio(u,v,w) = 0.0;
      end
    end
  end
end
name = strcat('SSRB_ratio_EmSub_noPost',filenameT,'.raw');
fid = fopen(name,'w');
fwrite(fid,SinoRatio,'float32');
fclose(fid);

%% perform median filtering of sino ratio
%for i=1:NSliTot
%  SinoRatio(:,:,i) = medfilt2(SinoRatio(:,:,i),[3 3]);
%end

%% Reconstruction of RatioComplete
%filenameMask  = 'SSRB_997ratio_mask_600_04TransPhantom2HOT.raw';
%filenameMask  = 'SSRB_mask_02TransPhantom2.raw';
%filenameMask  = 'SSRB_mask_03TransPhantom1.raw';
filenameMask  = 'mask_transPig1.raw';
%filenameMask  = 'mask_transPig2.raw';
SinoRatio = FillSinoGaps(name,filenameMask);
SinoRatio = smooth3(SinoRatio,'gaussian',[3 3 3],0.42466);
%SinoRatio = smooth3(SinoRatio,'gaussian',[9 9 7],2);

name = strcat('SSRB_ratio_SEG_',filenameT,'.raw');
fid = fopen(name,'w');
fwrite(fid,SinoRatio,'float32');
fclose(fid);

% rectangular masking of ratio sinogram
for i=1:Nbins
  for j=1:Nproj
    for k=1:NSliTot
      if i<122 || i>222 || k<25
        SinoRatio(i,j,k) = 0.0;
      end
    end
  end
end

recon = 4.83559*OSEM_Recon(SinoRatio,filenameT);
name = strcat('OSEM_ratio_SEG_',filenameT,'.raw');
fid  = fopen(name, 'w');
fwrite(fid,recon,'float32');
fclose(fid);

end