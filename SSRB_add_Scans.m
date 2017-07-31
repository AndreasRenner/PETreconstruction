function SSRB_add_Scans(filename)

% Basic Parameters of Siemens Biograph mMR
Nbins   = 344;         % Number of radial bins (NRAD)
Nproj   = 252;         % Number of projections (NANG)
Nslices = 127;         % Number of slices (2*Nrings -1)

% Number of sinogram parts within one scan
Nparts = 40;

% Add time from filling of Pellet to start of total scan
if     strcmp(filename, '05BlankFast')
  decayCor = 3471;
  Nscans = 10;
elseif strcmp(filename, '02TransPhantom1')
  decayCor = 2234;
  Nscans = 2;
elseif strcmp(filename, '03TransPhantom2')
  decayCor = 2630;
  Nscans = 2;
elseif strcmp(filename, '04TransRQ')
  decayCor = 3039;
  Nscans = 2;
end

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

% Load all Files into a list
name = strcat('sino_SSRB_',filename,'_*_*.raw');
d = dir(name);
filenames = {d.name};

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
  timespan(2*i-1) = difference(i)/2;
  timespan(2*i)   = difference(i)/2;  
end
for i=1:(Nparts*Nscans)
  meantime = 2000;
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
  
  for j=1:Nscans
    k = indices(j);
    totalCorF(k) = scanNumCor*timeF(k)/decayF(k);
    fprintf('Total Cor. Faktor for part %u, %u is: %u\r',i,j,totalCorF(k));
    name = cell2mat(filenames(k));
    fid  = fopen(name,'r');
    for l=1:Nslices
      Sino2D = fread(fid,[Nbins,Nproj],'float32');
      SIN2D  = double(Sino2D)*totalCorF(k);
      SSRBSino(:,:,l) = SSRBSino(:,:,l) + SIN2D;
    end
    fclose(fid);
  end
  
  % (4) SENSITIVITY CORRECTION
  % calculate total number of cor. counts
  totalcounts = sum(sum(sum(SSRBSino)));
  % apply correction factor
  sensitivityCor = 200000000./totalcounts;
  fprintf('Total Counts for part %u is: %u; corF: %u\r',i,totalcounts,sensitivityCor);
  SSRBSino = SSRBSino*sensitivityCor;
  
  name = strcat('SSRB_cor_',filename,'_',num2str(i),'.raw');
  fid = fopen(name,'w');
  fwrite(fid,SSRBSino,'float32');
  fclose(fid);
  
  for u=1:Nslices
    recon(:,:,u)=iradon(SSRBSino(:,:,u),theta,Nbins);
  end

  name = strcat('recon_', filename,'_',num2str(i),'.raw');
  fid = fopen(name, 'w');
  fwrite(fid,recon,'float32');
  fclose(fid);
end

end

