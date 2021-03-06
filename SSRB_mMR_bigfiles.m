function SSRB_mMR_bigfiles(filename,fileending,Nscans,ScatterCor)

% Estimate Number of Tags from Filesize
file  = strcat(filename, fileending);
size  = dir(file);
Ncs   = ceil(size.bytes/4);
fprintf('%s - Total Number of Tags: \t %10.0f\r',filename,Ncs);

Npos  = Ncs/Nscans;
%Nread = Npos;
%Nread = ceil(Npos*0.93);
%Nread = ceil(Npos*1.1);  % for transPig2
%Nread = ceil(Npos*1.11); % for blankPig2
%Nread = ceil(Npos*1.2);  % for transPig1
%Nread = ceil(Npos*1.14); % for blankPig1
%Nread = ceil(Npos*1.2);  % for transHP
Nread = ceil(Npos*1.14); % for blankHP
%Nread = ceil(Npos*1.15); % for 04HOT
%Nread = ceil(Npos*1.05); % for Trans
%Nread = ceil(Npos*1.16); % for 01Blank 
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
%scanDirection=1;
%
% Values for 03TransPhantom1 (2017-08-11)
%start = 100606+108831+106818+101711+104281+100627;
%tstart=[20580,123140,225080,326680,428970,534140,636800];
%tstop =[99290,202500,303510,406050,507720,613580,716540];
%faktor=[1.05,1.02,1,0.98,0.96,1,1]; %0 = 1.05
%scanDirection=1;
%
% Values for 04TransPhantom2HOT (2017-08-11)
%start = 106048+112690;%+279505+279352+111127+101297;
%tstart=[14440,120060,221710,679020,782860,889820,991260];
%tstop =[94090,198470,300740,757850,862400,968080,1070510];
%faktor=[0.9,1.5,1.42,0.82,0.75,1,1]; %0 = 0.9
%scanDirection=0;
% For Emission-Scan:
%start = 106048+112690+93222; % 317977
%faktor= 2*0.9 + 0.82; %0 = 1.15;
%
% Values for 01Blank (2017-08-11)
%start = 125764+89206+223582+99552+113780+83606+124739+91346+ ...
%    141472+113140+113379+86136+108868+103166;%+113781;
%tstart=[27730,133550,333310,441450,549910,652790,761150,869510, ...
%    984800,1094970,1206910,1323600,1426750,1532420,1635200];
%tstop=[107900,213860,412540,521480,629770,734250,841870,950490, ...
%   1063860,1174680,1286040,1403280,1505830,1611940,1712770];
%faktor=[1.14,1.11,1.08,1.07,1.05,1.06,1.01,1.02,1,0.94,0.92,0.92,0.9,0.88,1];
%scanDirection=1;
%
% Values for blankPig2 (2018-02-07)
%start = 79735+102648+82798+69419+99641+82791+80860+68081+95917+ ...
%    76789+81332+77737+69440+82514+86093+68218+85039+77650+69030+ ...
%    78518+93441+92075;
%tstart=[12990,102230,182650,265660,350020,434680,517280,598960, ...
%    679510,762400,839430,920980, 998410,1082450,1159190,1236640, ...
%    1313250,1390520,1467540,1545830,1624750,1708440];
%tstop =[78610,168480,248270,332140,415790,500590,582610,665590, ...
%    745580,829000,905500,987680,1066480,1149110,1225910,1303460, ...
%    1379680,1457210,1534920,1613630,1690610,1775880];
%faktor=[1.12,1.1,1.07,1.08,1.04,1.03,1.03,1.04,1,1,0.99,0.99,0.98, ...
%    0.97,0.94,0.94,0.93,0.92,0.93,0.9,1,1];
%scanDirection=0;
%
% Values for transPig2 (2018-02-07)
%start = 137238+167157+114013+156948+139553+138200+137543+140376+ ...
%    148679+133161+140618+246509;
%tstart=[ 70090,211610,350830,492190,630130,772870,911860,1052350, ...
%    1191050,1334250,1474950,1615950];
%tstop =[136270,278160,417420,558950,696540,839900,978430,1119220, ...
%    1256660,1401150,1540690,1682500];
%faktor=[1.11,1.04,1.06,1.02,1,0.98,0.96,0.94,0.92,0.9,1,1];
%scanDirection=0;
%
% Values for blankPig1 (2018-02-06)
%start = 75468+90774+104388+68421+106921+67879+101581+85532+86058+ ...
%    78595+81316+67406+93015+81387+80486+77879+67348+89667+77146+ ...
%    85348+92752;
%tstart=[ 8920, 99760,184280,271330,354880,446860,530660,616280, ...
%    702450,787640,866290, 947420,1027420,1107580,1188850,1269330, ...
%    1347920,1430910,1514810,1598750,1676270,incomplete];
%tstop =[74390,165680,248960,337370,419060,513590,594880,684160, ...
%    766600,854060,931980,1014250,1091610,1173860,1252880,1336140, ...
%    1413470,1497660,1580330,1664910,1741250,incomplete];
%faktor=[1.14,1.12,1.09,1.08,1.08,1.06,1.08,1.03,1.03,1.03,1,1,0.98, ...
%    0.96,0.96,0.94,0.96,0.92,0.93,0.9,1,1];
%scanDirection=0;
%
% Values for transPig1 (2018-02-06)
%start = 164280+127646+217956+103474+169664+115390+124558;%+185704+ ...
%    160555+70173+167843+192728;
%tstart=[ 67920,222670,369760,519820,659710,802400, 947150,1091270, ...
%    1229760,1371270,1511110,1655120];
%tstop =[138520,290160,435040,584770,723360,869840,1015760,1157110, ...
%    1294040,1437890,1576290,1722800];
%faktor=[1.1,1.1,1,1,1,1,0.98,0.92,0.88,0.9,1,1];
%scanDirection=0;
%
% Values for blankHP (2018-02-07)
%start = ;
%tstart=[];
%tstop =[];
%faktor=[];
%scanDirection=0;
%
% Values for transHP (2018-02-07)
%start = ;
%tstart=[];
%tstop =[];
%faktor=[];
%scanDirection=0;
%

fid   = fopen(file,'r');

%offset= uint64(2*ceil(Npos*1.14)+ceil(Npos*1.12)+ceil(Npos*1.09)+ ...
%    3*ceil(Npos*1.08)+ceil(Npos*1.06)+3*ceil(Npos*1.03)+2*Npos+ ...
%    ceil(Npos*0.98)+3*ceil(Npos*0.96)+ceil(Npos*0.94)+ ...
%    ceil(Npos*0.92))*4;%+ ...
%    ceil(Npos*1.07)+ceil(Npos*1.08)+ceil(Npos*1.04)+2*ceil(Npos*1.03)...
%    +ceil(Npos*1.04)+2*Npos+2*ceil(Npos*0.99)+ceil(Npos*0.98)+ ...
%    ceil(Npos*0.97)+2*ceil(Npos*0.94)+ceil(Npos*0.93)+ceil(Npos*0.92)+ ...
%    ceil(Npos*0.93))*4; %...
%fseek(fid,offset,'bof');

for i=1:Nscans
  dlist = fread(fid,[Nread],'uint32');
  
  % Get acquisition time in ms and output frequency of Tags
  fprintf('\n');
  fprintf('Information about Part %u:\r',i);
  acqTime = tagFrequency(dlist);
  fprintf('Acquisition Time Part %u: %u [ms]\r',i , acqTime);
  % Output prompts per second/10
  temporalPromptDistribution = p_s(dlist,acqTime,start);
  figure('Name', 'Temporal Distribution of Prompts');
  plot(temporalPromptDistribution);
  xlabel('Time [s/100]');
  ylabel('Prompts per second');
  
  if ScatterCor
    %fprintf('Input Scan Direction\r');
    %prompt='1 -> highest peak first; 0 -> else: ';
    %scanDirection=input(prompt);
    %clear prompt;
    % 01Blank         starts with highest peak first
    % 02TransPhantom2 starts with highest peak first
    % 03TransPhantom1 starts with highest peak first
    % 04TransPhantom2HOT starts with low  peak first
    % blankPig2          starts with low  peak first
    % transPig2          starts with low  peak first
    % transPig1          starts with low  peak first
    
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
    cutlist=zeros(1,161);
    overlap=zeros(1,161);
    cutlist(161)=minlist(21);
    for j=1:20
      step = (minlist(j+1)-minlist(j))/8;
      cutlist(8*j-7) = minlist(j);
      cutlist(8*j-6) = minlist(j)+step;
      cutlist(8*j-5) = minlist(j)+2*step;
      cutlist(8*j-4) = minlist(j)+3*step;
      cutlist(8*j-3) = minlist(j)+4*step;
      cutlist(8*j-2) = minlist(j)+5*step;
      cutlist(8*j-1) = minlist(j)+6*step;
      cutlist(8*j)   = minlist(j)+7*step;
      overlap(8*j-7) = step*0.1;
      overlap(8*j-6) = step*0.1;
      overlap(8*j-5) = step*0.1;
      overlap(8*j-4) = step*0.1;
      overlap(8*j-3) = step*0.1;
      overlap(8*j-2) = step*0.1;
      overlap(8*j-1) = step*0.1;
      overlap(8*j)   = step*0.1;
    end
    for j=1:160
      cut1 = cutlist(j)-overlap(j);
      cut2 = cutlist(j+1)+overlap(j+1);
      dlistcut = cutlmdata(dlist,cut1,cut2);
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
        dlistname = strcat(filename,stamp,num2str(161-j));
      else
        dlistname = strcat(filename,stamp,num2str(j));
      end
      sino = makeSinoRandomSubstraction(dlistcut);
      clear dlistcut;
      SSRB_mMR(sino,dlistname);
    end
    
    if scanDirection
      scanDirection = 0;
    else
      scanDirection = 1;
    end
    
  else
    % Get acquisition time in ms and output frequency of Tags
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
    %prompt = 'Enter start time in ms:';
    %tstart = input(prompt);
    %clear prompt;
    %prompt = 'Enter stop time in ms:';
    %tstop  = input(prompt);
    %clear prompt;

%    dlistcut  = cutlmdata(dlist,tstart(i),tstop(i));
%    clear dlist;
%
%    name = strcat(filename, num2str(i),'randomSubstracted');
%    % make Sino without random substraction
%    %sino = makeSinoBig(dlistcut);
%    % make Sino with random substraction
%    sino = makeSinoRandomSubstraction(dlistcut);
%    clear dlistcut;
%    SSRB_mMR(sino,name);
  end
  
  prompt = 'Enter Nread Faktor:';
  faktor = input(prompt);
  %Nread  = Npos;
  Nread  = ceil(Npos*faktor);
  %Nread  = ceil(Npos*faktor(i));
  
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
