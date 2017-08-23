function SSRB_mMR_bigfiles(filename, fileending, Nscans, ScatterCor)

% Estimate Number of Tags from Filesize
file  = strcat(filename, fileending);
size  = dir(file);
Ncs   = ceil(size.bytes/4);
fprintf('%s - Total Number of Tags: \t %10.0f\r',filename,Ncs);

Npos  = Ncs/Nscans;
%Nread = Npos;
Nread = ceil(Npos*1.16); % for 01Blankceil(Npos*1.05); % for Trans 
fprintf('Number of Tags per read: %u\r\n', Nread);
start = 0;
%
% Values for 05BlankFast (2017-05-19)
%start = 154938+118585+110162+155863+87207+124791; %97032+157230+96259
%tstart=[55380,177360,297590,422600,540230,658250,768720,879910,1006230,1115500];
%tstop=[140530,259620,382720,506210,625560,740220,848090,971760,1087890,1197180];
%faktor=[1.05,1.05,1.05,1,0.98,0.9,1.08,0.9,1,1];
%
% Values for 02TransPhantom2 (2017-08-11)
%start = 122850+110091+100678+104987+118226;
%tstart=[27440,131270,236000,336690,440090,585420];
%tstop=[106750,211540,314470,416200,518670,664460];
%faktor=[1.03,0.97,0.98,0.95,0.95,1]; %0 = 1.05
%
% Values for 03TransPhantom1 (2017-08-11)
%start = 100606+108831+106818+101711+104281+100627;
%tstart=[20580,123140,225080,326680,428970,534140,636800];
%tstop =[99290,202500,303510,406050,507720,613580,716540];
%faktor=[1.05,1.02,1,0.98,0.96,1,1]; %0 = 1.05
%
% Values for 01Blank (2017-08-11)
%start = 125764+89206+223582+99552+113780+83606+124739+91346+ ...
%    141472+113140+113379+86136+108868+103166;%+113781;
%tstart=[27730,133550,333310,441450,549910,652790,761150,869510, ...
%    984800,1094970,1206910,1323600,1426750,1532420,1635200];
%tstop=[107900,213860,412540,521480,629770,734250,841870,950490, ...
%   1063860,1174680,1286040,1403280,1505830,1611940,1712770];
%faktor=[1.14,1.11,1.08,1.07,1.05,1.06,1.01,1.02,1,0.94,0.92,0.92,0.9,0.88,1];
fid   = fopen(file,'r');
%offset= uint64(ceil(Npos*1.16)+ceil(Npos*1.14)+ceil(Npos*1.11)+ ...
%               ceil(Npos*1.08)+ceil(Npos*1.07)+ceil(Npos*1.05)+ ...
%               ceil(Npos*1.06)+ceil(Npos*1.01)+ceil(Npos*1.02)+ ...
%               Npos+ceil(Npos*0.94)+2*ceil(Npos*092)+ ...
%               ceil(Npos*0.9))*4;%+ceil(Npos*0.88))*4;
%fseek(fid,offset,'bof');

for i=1:(Nscans-1)
  dlist = fread(fid,[Nread],'uint32');
  
  % Get acquisition time in ms and output frequency of Tags
  fprintf('\n');
  fprintf('Information about Part %u:\r',i);
  acqTime = tagFrequency(dlist);
  fprintf('Acquisition Time Part %u: %u [ms]\r',i , acqTime);
  % Output prompts per second/10
  %temporalPromptDistribution = p_s(dlist,acqTime,start);
  %figure('Name', 'Temporal Distribution of Prompts');
  %plot(temporalPromptDistribution);
  %xlabel('Time [s/100]');
  %ylabel('Prompts per second');
  
  if ScatterCor
    fprintf('Input Scan Direction\r');
    prompt='1 -> highest peak first; 0 -> else: ';
    scanDirection=input(prompt);
    clear prompt;
    name = strcat('minlistScan',num2str(i),filename);
    try
      load(name);
    catch ME
      if (strcmp(ME.identifier,'MATLAB:load:couldNotReadFile'))
        warning('Did not find >>minlistBlanc.mat<<. Enter Minima manually.');
        % Get acquisition time in ms and output frequency of Tags
        minlist=zeros(1,21);
        for j=1:21
          fprintf('%u. Minimum - ',j);
          prompt='Enter time in ms:';
          minlist(j)=input(prompt);
        end
        clear prompt;
        save(name,'minlist');
      else
        rethrow(ME)
      end
    end
    cutlist=zeros(1,41);
    cutlist(41)=minlist(21);
    for j=1:20
      temp = (minlist(j+1)-minlist(j))/2;
      cutlist(2*j-1) = minlist(j);
      cutlist(2*j)   = minlist(j)+temp;
    end
    for j=1:40
      dlistcut = cutlmdata(dlist,cutlist(j),cutlist(j+1));
      % Check if Cut was succesful
      %acqTime  = tagFrequency(dlistcut);
      %tempDist = p_s(dlistcut,acqTime,cutlist(j));
      %figure('Name', 'Prompt Distribution of Cutlist');
      %plot(tempDist);
      %xlabel('Time [s/100]');
      %ylabel('Prompts per second');
      timestamp = cutlist(j)+(cutlist(j+1)-cutlist(j))/2;
      stamp = strcat('_',num2str(timestamp),'_');
      if scanDirection
        dlistname = strcat(filename,stamp,num2str(41-j));
      else
        dlistname = strcat(filename,stamp,num2str(j));
      end
      sino = makeSinoBig(dlistcut);
      clear dlistcut;
      SSRB_mMR(sino,dlistname);
    end
    
  else
    %prompt = 'Enter start time in ms:';
    %tstart = input(prompt);
    %clear prompt;
    %prompt = 'Enter stop time in ms:';
    %tstop  = input(prompt);
    %clear prompt;

    dlistcut  = cutlmdata(dlist,tstart(i),tstop(i));
    clear dlist;

    name = strcat(filename, num2str(i),'randomSubstracted');
    % make Sino without random substraction
    %sino = makeSinoBig(dlistcut);
    % make Sino with random substraction
    sino = makeSinoRandomSubstraction(dlistcut);
    clear dlistcut;
    SSRB_mMR(sino,name);
  end
  
  %prompt = 'Enter Nread Faktor:';
  %faktor = input(prompt);
  %Nread  = Npos;
  Nread  = ceil(Npos*faktor(i));
  
  start  = start+acqTime;
end
fclose(fid);

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
% Make uncompressed Sinogram without random substraction
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
