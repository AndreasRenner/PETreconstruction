function SSRB_mMR_emission(filename,fileending,Nscans)

% Estimate Number of Tags from Filesize
file  = strcat(filename, fileending);
size  = dir(file);
Ncs   = ceil(size.bytes/4);
fprintf('%s - Total Number of Tags: \t %10.0f\r',filename,Ncs);

Npos  = Ncs/Nscans;
Nread = Npos;
%Nread = ceil(Npos*0.93);
%Nread = ceil(Npos*1.1);  % for transPig2
%Nread = ceil(Npos*1.2);  % for transPig1
%Nread = ceil(Npos*1.2);  % for transHP
%Nread = ceil(Npos*1.15); % for 04HOT
fprintf('Number of Tags per read: %u\r\n', Nread);
%start = 0;
%
% Values for 04TransPhantom2HOT (2017-08-11)
%tstart=[14440,120060,221710,679020,782860,889820,991260];
%tstop =[94090,198470,300740,757850,862400,968080,1070510];
%faktor=[0.9,1.5,1.42,0.82,0.75,1,1]; %0 = 0.9
% For Emission-Scan:
%faktor= 2*0.9 + 0.82; %0 = 1.15;
%
% Values for transPig2 (2018-02-07)
%tstart=[ 0,141,283,423,563,703,844, 982,1123,1263,1406,1546,1688];
%tstop =[64,205,347,487,627,767,908,1046,1187,1328,1470,1610,1799];
%faktor=[1,1,1,1,1,1,1,1,1,1,3,1];
%
% Values for transPig1 (2018-02-06)
tstart=[ 0,146,299,444,591,733,876,1021,1162,1302,1443,1586,1737];
tstop =[63,209,362,507,654,796,940,1084,1225,1365,1507,1649,1799];
faktor=[1,1,1,1,1,1,1,1,1,1,1.4,1];
%
% Values for transHP (2018-02-07)
%tstart=[];
%tstop =[];
%faktor=[];

fid   = fopen(file,'r');

for i=1:Nscans
  dlist = fread(fid,[Nread],'uint32');
  
%  % Get acquisition time in ms and output frequency of Tags
%  fprintf('\n');
%  fprintf('Information about Part %u:\r',i);
%  acqTime = tagFrequency(dlist);
%  fprintf('Acquisition Time Part %u: %u [ms]\r',i , acqTime);
%  % Output prompts per second/10
%  temporalPromptDistribution = p_s(dlist,acqTime,start);
%  figure('Name', 'Temporal Distribution of Prompts');
%  plot(temporalPromptDistribution);
%  xlabel('Time [s/100]');
%  ylabel('Prompts per second');

  dlistcut = cutlmdata(dlist,tstart(i)*1000,tstop(i)*1000);
  if i<Nscans
    clear dlist;
  end
  name = strcat(filename, num2str(i),'emission');
  % make Sino without random substraction
  %sino = makeSinoBig(dlistcut);
  % make Sino with random substraction
  sino = makeSinoRandomSubstraction(dlistcut);
  clear dlistcut;
  SSRB_mMR(sino,name);
  
  %prompt = 'Enter Nread Faktor:';
  %faktor = input(prompt);
  %Nread  = ceil(Npos*faktor);
  Nread = ceil(Npos*faktor(i));
%  start  = start+acqTime;
end
fclose(fid);

dlistcut = cutlmdata(dlist,tstart(i+1)*1000,tstop(i+1)*1000);
clear dlist;
name = strcat(filename, num2str(i+1),'emission');
sino = makeSinoRandomSubstraction(dlistcut);
clear dlistcut;
SSRB_mMR(sino,name);

end

% -------------------------------------------------------
% Make uncompressed Sinogram without random substraction
function sino = makeSinoBig(dlist)
% Basic Parameters of Siemens Biograph mMR
Nbins   = 344;         % Number of radial bins (NRAD)
Nproj   = 252;         % Number of projections (NANG)
Nplanes = 4084;        % Number of 3D sinogram planes (Nsinos)
sinoDim = Nbins * Nproj * Nplanes; % 354 033 792

% Prompts between [001111...1] and [01111...1]
ptag = dlist((dlist<2^31)&(dlist>=2^30));
ptag = ptag-(2^30);
% clear out-of-sinograms indices
% sinoDim equal to 1 0101 0001 1010 0010 0000 1000 0000
ptag = ptag((ptag<=prod(sinoDim))&(ptag>0));
% build sinograms
sinolist = accumarray(ptag,1,[sinoDim,1]);
clear ptag;
index = 0;
for k=1:Nplanes
  for j=1:Nproj
    for i=1:Nbins
      index = index+1;
      sino(i,j,k) = sinolist(index);
    end
  end
end
end

% -------------------------------------------------------
% Make uncompressed Sinogram with random substraction
function sino = makeSinoRandomSubstraction(dlist)
% Basic Parameters of Siemens Biograph mMR
Nbins   = 344;         % Number of radial bins (NRAD)
Nproj   = 252;         % Number of projections (NANG)
Nplanes = 4084;        % Number of 3D sinogram planes (Nsinos)
sinoDim = Nbins * Nproj * Nplanes; % 354 033 792

% Prompts between [001111...1] and [01111...1]
ptag = dlist((dlist<2^31)&(dlist>=2^30));
ptag = ptag-(2^30);
% clear out-of-sinograms indices
% sinoDim equal to 1 0101 0001 1010 0010 0000 1000 0000
ptag = ptag((ptag<=prod(sinoDim))&(ptag>0));
% build sinograms
sinolist = accumarray(ptag,1,[sinoDim,1]);
clear ptag;

% Randoms lower or equal [001111...1]
rtag = dlist((dlist<2^30));
% clear out-of-sinograms indices
rtag = rtag((rtag<=prod(sinoDim))&(rtag>0));
% substract sinograms
sinolist = sinolist-accumarray(rtag,1,[sinoDim,1]);
clear rtag;

index = 0;
for k=1:Nplanes
  for j=1:Nproj
    for i=1:Nbins
      index = index+1;
      sino(i,j,k) = sinolist(index);
    end
  end
end
end


% -------------------------------------------------------
% Cut the first X and after Y ms of acquisition
function dlist = cutlmdata(dlist,tstart,tstop)
% Convert time in [ms] to the format of a Ttag and search for
% the corresponding position of the Ttag in dlist 
posx = find(dlist==(int64(tstart)+2^31));
posy = find(dlist==(int64(tstop)+2^31));
  
% Cut data list
dlist = dlist(posx:posy);  
end
