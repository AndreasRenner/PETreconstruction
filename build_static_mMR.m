function build_static_mMR(filenameBlanc, filenameTrans)
%filenameBlanc='BlancScan.IMA';
%filenameTrans='TransmissionScan.IMA';

% BLANC-SCAN: read file and make sinogram
dlistB = readfile(filenameBlanc);
% Load minima from file or querry user to enter minima
try
  load('minlistBlanc.mat');
catch ME
  if (strcmp(ME.identifier,'MATLAB:load:couldNotReadFile'))
    warning('Did not find >>minlistBlanc.mat<<. Enter Minima manually.');
    % Get acquisition time in ms and output frequency of Tags
    acqTimeB = tagFrequency(dlistB);
    fprintf('BlancScan - Acquisition time: %u [ms]\r', acqTimeB);
    % Output prompts per second/10
    % ->Output should be used to determine the time when
    %   Transmission source is entering/leaving helical path
    acqTimeSecB = floor(acqTimeB*0.01);
    temporalPromptDistributionB = p_s(dlistB,acqTimeSecB);
    %findpeaks(temporalPromptDistribution,'Annotate','extents','WidthReference','halfheight');
    figure('Name','BlancScan: Temporal Distribution of Prompts');
    plot(temporalPromptDistributionB);
    xlabel('Time [s/10]');
    ylabel('Prompts per second');
    minlistBlanc=zeros(1,21);
    for i=1:21
      fprintf('%u. Minimum - ',i);
      prompt='Enter time in ms:';
      minlistBlanc(i)=input(prompt);
    end
    clear prompt;
  else
    rethrow(ME)
  end
end
[dlistB] = cutlmdata(dlistB,minlistBlanc(1),minlistBlanc(21));
makeSino(dlistB,filenameBlanc);
clear dlistB;
SSRB_Blanc =  SSRB_mMR(filenameBlanc);

% TRANSMISSION-SCAN: read file and make sinogram
dlistT = readfile(filenameTrans);
% Load minima from file or querry user to enter minima
try
  load('minlistTrans.mat');
catch ME
  if (strcmp(ME.identifier,'MATLAB:load:couldNotReadFile'))
    warning('Did not find >>minlistTrans.mat<<. Enter Minima manually.');
    % Get acquisition time in ms and output frequency of Tags
    acqTimeT = tagFrequency(dlistT);
    fprintf('TransScan - Acquisition time: %u [ms]\r', acqTimeT);
    % Output prompts per second/10
    acqTimeSecT = floor(acqTimeT*0.01);
    temporalPromptDistributionT = p_s(dlistT,acqTimeSecT);
    figure('Name','TransScan: Temporal Distribution of Prompts');
    plot(temporalPromptDistributionT);
    xlabel('Time [s/10]');
    ylabel('Prompts per second');
    minlistTrans=zeros(1,21);
    for i=1:21
      fprintf('%u. Minimum - ',i);
      prompt='Enter time in ms:';
      minlistTrans(i)=input(prompt);
    end
    clear prompt;
  else
    rethrow(ME)
  end
end
[dlistT] = cutlmdata(dlistT,minlistTrans(1),minlistTrans(21));
makeSino(dlistT,filenameTrans);
clear dlistT;
SSRB_Trans = SSRB_mMR(filenameTrans);

% Calculate and reconstruc ratio of Blank- to Transmission-Scan
reconSinoRatio(SSRB_Blanc,SSRB_Trans);
end

% -------------------------------------------------------
% Read file into dlist
function dlist = readfile(filename)
  % Estimate Number of Tags from Filesize
  size = dir(filename);
  Ncs  = ceil(size.bytes/4);
  fprintf('%s - Number of Tags: \t %10.0f\r\n',filename,Ncs);
  
  fid  = fopen(filename,'r');
  dlist= fread(fid,[Ncs],'uint32');    
  fclose(fid);
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
