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
SSRB_Blanc = makeSino(dlistB,filenameBlanc);
clear dlistB;

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
SSRB_Trans = makeSino(dlistT,filenameTrans);
clear dlistT;

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
% Get frequency of all individual Tags
function duration = tagFrequency(dlist)
  % Event Packets:  0XXX ...
  prompts    = 0; % 01XX ...
  delays     = 0; % 00XX ...
  % Tag Packets:    1XXX ...
  timeTag    = 0; % 100X ...
  DtimeTag   = 0; % 101X ...
  motionTag  = 0; % 1100 ...
  patientTag = 0; % 1110 ...
  controlTag = 0; % 1111 ...
  currentTime= 0;
  
  fprintf('****************************************************************\r');
  fprintf('START Detailed Information about Frequency of individual Tags:\r');
  fprintf('Special Control/Patient Tags:\r');
  specialTagCounter = 0;
  
  for i=1:length(dlist)
    if       dlist(i) < 1073741824
      delays = delays + 1;
      %fprintf('Pos %u\tat %s ms:\tDelay Event:\t%s\r',i,currentTime,num2str(dec2bin(dlist(i))));
    elseif   dlist(i) < 2147483648
      prompts = prompts + 1;
      %fprintf('Pos %u\tat %s ms:\tPrompt Event:\t%s\r',i,currentTime,num2str(dec2bin(dlist(i))));
    elseif   dlist(i) < 2684354560
      timeTag = timeTag + 1;
      currentTime=num2str(dlist(i)-2^31);
      %fprintf('Pos %u\tat %s ms:\tTime Tag:\t%s\r',i,currentTime,num2str(dec2bin(dlist(i))));
    elseif   dlist(i) < 3221225472
      DtimeTag = DtimeTag + 1;
      %fprintf('Pos %u\tat %s ms:\tDead-Time Tag:\t%s\r',i,currentTime,num2str(dec2bin(dlist(i))));
    elseif   dlist(i) < 3758096384
      motionTag = motionTag + 1;
      %fprintf('Pos %u\tat %s ms:\Motion Tag:\t%s\r',i,currentTime,num2str(dec2bin(dlist(i))));
    elseif   dlist(i) < 4026531840
      patientTag = patientTag + 1;
      if     dlist(i)== 3774877696
        specialTagCounter = specialTagCounter + 1;
        fprintf('Pos %u\tat %s ms:\tPatientTag - Respiratory Trigger Gate on\r',i,currentTime);
      elseif dlist(i)== 3774877697
        specialTagCounter = specialTagCounter + 1;
        fprintf('Pos %u\tat %s ms:\tPatientTag - Respiratory Trigger Gate off\r',i,currentTime);
      end
      %fprintf('Pos %u\tat %s ms:\tPatient Tag:\t%s\r',i,currentTime,num2str(dec2bin(dlist(i))));
    else
      controlTag = controlTag + 1;
      if     dlist(i)==4294901760
        fprintf('Pos %u\tat %s ms:\tControlTag - Start of Acquisition\r',i,currentTime);
      elseif dlist(i)==4286545920
        fprintf('Pos %u\tat %s ms:\tControlTag - Time Synchronization with MR\r',i,currentTime);
      end
      %fprintf('Pos %u\tat %s ms:\tControl Tag:\t%s\r',i,currentTime,num2str(dec2bin(dlist(i))));     
    end
  end
  clear dlist;
  if specialTagCounter==0
    fprintf('None\r');
  end
  
  events = prompts + delays;
  tags   = timeTag + DtimeTag + motionTag + patientTag + controlTag;
  fprintf('Total number of Event Packets: \t %10.0f\r', events);
  fprintf('Total number of Tag Packets:   \t %10.0f\r', tags);
  fprintf('Prompt Events:  \t %10.0f\r', prompts);
  fprintf('Delay Events:   \t %10.0f\r', delays);
  fprintf('Time Tags:      \t %10.0f\r', timeTag);
  fprintf('Dead-Time Tags: \t %10.0f\r', DtimeTag);
  fprintf('Motion Tags:    \t %10.0f\r', motionTag);
  fprintf('Patient Tags:   \t %10.0f\r', patientTag);
  fprintf('Control Tags:   \t %10.0f\r', controlTag);
  fprintf('END Detailed Information about Frequency of individual Tags\r');
  fprintf('****************************************************************\r\n');
  
  % Return scan-duration in [ms]
  duration = timeTag;
end

% -------------------------------------------------------
% Get prompts per second/10
function [z] = p_s(dlist,time)
  p=0;
  T=0;
  tmplen=0;
  z = zeros(time,1);
  for i=1:length(dlist)
    if (dlist(i)<(2^31))&&(dlist(i)>=(2^30))
      p=p+1;
    elseif (dlist(i)>=(2^31))&&(dlist(i)<2684354560)
      Ttag=dlist(i);
      if mod(Ttag-(2^31),100)==0
        tmptime = T*0.01;
        if tmptime>0
          z(tmptime) = p-tmplen;
          tmplen=p;
        end
      end
      T=T+1;
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

% -------------------------------------------------------
% Create Sinograms and perform SSRB (Single-Slice-ReBining)
function [SSRBSino] = makeSino(dlist,filename)
  % Basic Parameters of Siemens Biograph mMR
  Nbins   = 344;         % Number of radial bins (NRAD)
  Nproj   = 252;         % Number of projections (NANG)
  Nplanes = 4084;        % Number of 3D sinogram planes (Nsinos)
  Nseg    = 121;         % Number of segments
  Nrings  = 64;          % Number of Detector-Rings
  Nslices = 127;         % Number of slices (2*Nrings -1)
  % The Sinogram is arranged in 121 Segments; the first Segment consists of
  % the 64 direct planes, 2nd and 3rd Segment of 63 planes each for +1 and
  % -1 rings and so on; last and second-last Segment consist of 4 planes
  % each as the maximum ring difference is 60; this sums up to a total of
  % 4084 planes.
  sinoDim = Nbins * Nproj * Nplanes;
  Span    = 1;

  % Prompts between [001111...1] and [01111...1]
  ptag = dlist((dlist<2^31)&(dlist>=2^30));
  ptag = ptag-(2^30);
  % clear out-of-sinograms indices
  if ptag(ptag>sinoDim)
    ptag = ptag((ptag<=sinoDim));
  elseif ptag(ptag<=0)
    ptag = ptag((ptag>0));
  end
  % build sinograms
  sino = accumarray(ptag,1,[sinoDim,1]);
  clear ptag;
  
  % Randoms lower or equal [001111...1]
  rtag = dlist((dlist<2^30));
  % clear out-of-sinograms indices
  if rtag(rtag>sinoDim)
    rtag = rtag((rtag<=sinoDim));
  elseif rtag(rtag<=0)
    rtag = rtag((rtag>0));
  end
  % substract sinograms
  sino = sino-accumarray(rtag,1,[sinoDim,1]);
  clear rtag;
  
  % write sinograms to file
  sinogramname = strcat('sino_', filename, '.raw');
  fid=fopen(sinogramname,'w');
  fwrite(fid,uint16(sino),'uint16');
  fclose(fid);
  
  % SSRB Algorithm: creates a rebinned direct sinogram  
  SSRBSino = zeros(Nbins,Nproj,Nslices,'double');
  % INCLINATION of each segment
  Offset(1)=0;
  Offset(2)=(Span+1)/2;
  Offset(3)=(Span+1)/2;
  for iseg=4:Nseg
    Offset(iseg)=Offset(iseg-2)+Span;
  end
  % Number of sinograms in each segment
  Nsinos=0;
  for iseg=1:Nseg
    NsinoSeg(iseg)=Nrings-Offset(iseg);
    Nsinos=Nsinos+NsinoSeg(iseg);
  end

  fid=fopen(sinogramname,'r');
  Nsino=0;
  for iseg=1:Nseg
    for u=1:NsinoSeg(iseg)
      Nsino  = Nsino+1;
      Sino2D = fread(fid, [Nbins, Nproj], '*int16');
      k_SSRB = (2*u-1) + Offset(iseg);
      SIN2D  = double(Sino2D); 
      SSRBSino(:,:,k_SSRB)=SSRBSino(:,:,k_SSRB)+SIN2D; 
    end
  end
  fclose(fid);
  clear Sino2D SIN2D;
  
  SINR = double(SSRBSino);
  
  % Gaps Finder and Gap Filling (using inpaint)
  MASK = zeros(Nbins,Nproj);
  for u=1:Nslices
    MASK = MASK + SINR(:,:,u); 
  end  
  ngap=0;
  nnogap=0;
  for ith=1:Nproj
    for ir=1:Nbins
      if (MASK(ir,ith)<0.01)     
        ngap=ngap+1;    
        MASK(ir,ith)=0.0;
      else
        nnogap=nnogap+1; 
        MASK(ir,ith)=1.0;
      end
    end
  end  
  for u=1:Nslices
    SIN=SINR(:,:,u);
    SIN(MASK==0.)=nan;
    SIN2=inpaint_nans(SIN,2);
    SSRBSino(:,:,u)=SIN2;
  end
  clear SINR;
  
  % Smooth 3D data with gaussian kernel
  %SSRBSino = smooth3(SSRBSino,'gaussian',[3 3 3],0.42466);
  % sd of 0.42466 correlates to FWHM of 1
  SSRBSino = smooth3(SSRBSino,'gaussian',[5 5 5],1.5);
  
  newName = strcat('sino_SSRB_', filename, '.raw');
  fid = fopen(newName,'w');
  fwrite(fid,SSRBSino,'float32');
  fclose(fid);

%   % Transform sinogram do STIR projectiondata
%   for u=1:Nslices
%     for v=1:Nproj
%       for w=1:Nbins
%         STIR_proj(w,u,v)=SSRBSino(w,v,u);
%       end
%     end
%   end
%   stirName = strcat('STIR_', filename, '.s');
%   fid = fopen(stirName, 'w');
%   fwrite(fid,STIR_proj,'float32');
%   fclose(fid);
end

% -------------------------------------------------------
% Reconstruct ratio between Blanc- and Transmission-Scan
function reconSinoRatio(SSRB_Blanc,SSRB_Trans)
  %load('SSRB_Blanc.mat');
  %load('SSRB_Trans.mat');
  %SSRB_Blanc=smooth3(SSRB_Blanc,'gaussian',[7 7 5],2);
  %SSRB_Trans=smooth3(SSRB_Trans,'gaussian',[7 7 5],2);
  for u=1:size(SSRB_Blanc,1)
    for v=1:size(SSRB_Blanc,2)
      for w=1:size(SSRB_Blanc,3)
        if SSRB_Trans(u,v,w)<0.1 % SSRB_Blanc(u,v,w)<4 &&
          SSRB_Ratio(u,v,w)=log(SSRB_Blanc(u,v,w));
        else
          SSRB_Ratio(u,v,w)=log(SSRB_Blanc(u,v,w)./SSRB_Trans(u,v,w));
        end
      end
    end
  end
  
  newName = 'sino_SSRB_ratio.raw';
  fid = fopen(newName,'w');
  fwrite(fid,SSRB_Ratio,'float32');
  fclose(fid);
  
  theta = 180./size(SSRB_Blanc,2);
  for i=1:size(SSRB_Blanc,3)
    recon(:,:,i)=iradon(SSRB_Ratio(:,:,i),theta);
  end
  
  fid = fopen('reconRatio.raw', 'w');
  fwrite(fid,recon,'float32');
  fclose(fid);
end
