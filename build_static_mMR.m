function build_static_mMR(filename)

% Choose the file from an UI
%[filename, path]=uigetfile('*.*');
%cd(path);
% Get the Number of Counts (Ncs) from the size of the file
%cd('/home/andreas/data/PET_raw_data_20160603');

% Estimate Number of Tags from Filesize
filesize = dir(filename);
Ncs      = ceil(filesize.bytes/4);
clear filesize;
fprintf('Estimated Number of Tags: \t %10.0f\r\n', Ncs);

% Read File into a List
fid   = fopen(filename,'r');
dlist = fread(fid,[Ncs],'uint32');    
fclose(fid);
clear fid;
%cd('/home/andreas/code/PETreconstruction');

% Get acquisition time in ms and output frequency of Tags
acqTime = tagFrequency(dlist);
fprintf('Acquisition time: %u [ms]\r', acqTime);

% Output prompts per second/10
% ->Output should be used to determine the time when
%   Transmission source is entering/leaving helical path
acqTimeSec = floor(acqTime*0.01);
temporalPromptDistribution = p_s(dlist,acqTimeSec);
figure('Name','Temporal Distribution of Prompts');
plot(temporalPromptDistribution);
xlabel('Time [s/10]');
ylabel('Prompts per second');
%findpeaks(temporalPromptDistribution,'Annotate','extents','WidthReference','halfheight');

% Query the user to enter the minima or alternatively
% load minima from a file
promptinit = 'Is there an existing list of Minima? [Y=1]';
if input(promptinit)==1
  load('minlist.mat');
  %load('minlistBlanc.mat');
else
  minlist=zeros(1,21);
  for i=1:21
    fprintf('%u. Minimum - ',i);
    prompt='Enter time in ms:';
    minlist(i)=input(prompt);
  end
  clear prompt;
end
clear promptinit;

figure('Name','Extracted Round-Times from Minima of Prompts/Second');
plot(minlist);
xlabel('Round Number');
ylabel('Acquisition Time in [ms]');

rounds = (1:21);
fitMinList = fit(rounds',minlist','poly1');
evalFitMin = feval(fitMinList,rounds);

% Cut the first X and the last Y ms of acquisition
[dlist] = cutlmdata(dlist,evalFitMin(1),evalFitMin(21));

% Output detailed dead-time information and
% Get TimeTag of SingleBuckets
[singleBucketTimes,singleBuckets] = deadTimeInfo(dlist);

% Transform PET-acquisition-time into Transmission-scan-time
singleBucketTimes = singleBucketTimes - evalFitMin(1);

% Show Bucket-Single rates
showBucketSingles(singleBuckets,singleBucketTimes,evalFitMin);

% Create Sinograms
makeSino(cutdlist);

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
% Get Single-Buckets and output detailed dead-time information
function [singleBucket,blockSingles] = deadTimeInfo(dlist)
  D=0;
  bucketRounds = 1;
  LostEventsFirstLossyNode=0;
  LostEventsSecondLossyNode=0;
  millionEvents=0;
  for i=1:length(dlist)
    if (dlist(i)>=(2^31))&&(dlist(i)<2684354560)
      % Time Tag
      Ttag=num2str(dlist(i)-2^31);
    elseif (dlist(i)>=2684354560)&&(dlist(i)<3221225472)
      % Dead-Time Tag
      D=D+1;
      binaryTag=num2str(dec2bin(dlist(i)));
      typefield=bin2dec(binaryTag(4:6));
      if typefield==7
        loss=bin2dec(binaryTag(13:32));
        LostEventsFirstLossyNode=LostEventsFirstLossyNode+loss;
        millionEvents=millionEvents+1;
        %fprintf('1. Lossy Node: %u\tlost events at %s [ms]\r',loss,Ttag);
      elseif typefield==6
        loss=bin2dec(binaryTag(13:32));
        LostEventsSecondLossyNode=LostEventsSecondLossyNode+loss;
        millionEvents=millionEvents+1;
        %fprintf('2. Lossy Node: %u\tlost events at %s [ms]\r',loss,Ttag);
      else
        blocknum=bin2dec(binaryTag(4:13));
        singles=bin2dec(binaryTag(14:32));        
        blockSingles(bucketRounds,blocknum+1) = singles;
        %fprintf('Block: %u\tSingles: %u\t Time[ms]: %s\r',blocknum,singles,Ttag);
        if blocknum == 223
          singleBucket(bucketRounds) = str2num(Ttag);
          bucketRounds = bucketRounds + 1;
        end
      end
    end
  end
  totalLoss=LostEventsSecondLossyNode+LostEventsFirstLossyNode;
  totalEvents=millionEvents*1048575;
  % 1048575 corresponds to the 20 bit ,,lost event tally field''
  fprintf('Total Dead-time marks: %u\r',D);
  fprintf('Number of lost events inserted by the "1. Lossy Node": %u\r', LostEventsFirstLossyNode);
  fprintf('Number of lost events inserted by the "2. Lossy Node": %u\r', LostEventsSecondLossyNode);
  fprintf('Total number of lost event packets: \t %u\r', totalLoss);
  fprintf('Total number of initial events: \t~%u\r', totalEvents);
end

% -------------------------------------------------------
% Create Sinograms
function makeSino(dlist)
  % Basic Parameters of Siemens Biograph mMR
  Nbins   = 344;         % Number of radial bins
  Nproj   = 252;         % Number of projections
  Nplanes = 4084;        % Number of 3D sinogram planes
  sinoDim = Nbins * Nproj * Nplanes;
  
  % Prompts between [001111...1] and [01111...1]
  ptag = dlist(find((dlist<2^31)&(dlist>=2^30)));
  % Randoms lower or equal [001111...1]
  rtag = dlist(find((dlist<2^30)));

  % remove preceding tag if necessary
  ptag = ptag-(2^30);
  % clear out-of-sinograms indices
  ptag = ptag((ptag<=sinoDim)&(ptag>0));
  rtag = rtag((rtag<=sinoDim)&(rtag>0));    

  % build sinograms
  sino=accumarray(ptag,1,[sinoDim,1])-accumarray(rtag,1,[sinoDim,1]);
    
  % write sinograms to file
  sinogramname = strcat('sinogram_static_', filename, '.raw');
  fid=fopen(sinogramname,'w');
  fwrite(fid,uint16(sino),'uint16');
  fclose(fid);

  % Output for user
  fprintf('Finished reading list file\r\n');
  fprintf('Prompts\t\t:\t%s\r',num2str(length(ptag)));
  fprintf('Randoms\t\t:\t%s\r',num2str(length(rtag)));
end

% -------------------------------------------------------
% Cut the first X and after Y ms of acquisition
function cutdlist = cutlmdata(dlist,tstart,tstop)
  % Read values for X and Y
  %promptx='Enter time in ms you want first cut:';
  %prompty='Enter time in ms you want second cut:';
  %tstart=input(promptx);
  %tstop=input(prompty);
  
  % Convert time in [ms] to the format of a Ttag and search for
  % the corresponding position of the Ttag in dlist 
  posx = find(dlist==(int64(tstart)+2^31));
  posy = find(dlist==(int64(tstop)+2^31));
    
  % Create the cut data list
  cutdlist = dlist(posx:posy);
end

% -------------------------------------------------------
% Show single-bucket rates
function showBucketSingles(BlockSingles,timeList,minlist)
  timesteps = length(timeList);
  minlist = minlist-minlist(1);
  tstop = minlist(21);
  
  % Show summary plot of Bucket-Singles
  colorMax = floor(max(max(BlockSingles))/1000)*1000;
  figure();
  contourf(BlockSingles,10);
  colormap(hot);
  caxis([0 colorMax]);
  c = colorbar;
  c.Label.String = 'Bucket Singles Rate';
  ylabel('Ring Number');
  xlabel('Time [s]');
  
  % Gaussian-Fit________________________________________
  backward = 0;
  for k = 1:timesteps
    z = zeros(28,8);
    readPosition = 1;
    for i=1:8
      for j=1:28
        z(j,i) = BlockSingles(k,readPosition);
        readPosition = readPosition + 1;
      end
    end
    %figure('Name', 'Initial Distribution of Single-Rates');
    %imagesc(z);
    zmin = min(min(z));
    zmax = max(max(z));
    downshift=0;
    upshift=0;
    if (zmax-zmin)<1000
      fprintf('%u: Peak is too small for fitting!\r',k);
    else
      [izmax,jzmax] = find(z==zmax);
      if k==1 && jzmax(1)>6
        backward=1;
      end
      if izmax<7
        temp = zeros(28,8);
        for i=1:22
          temp((i+6),:)=z(i,:);
        end
        for i=23:1:28
          temp((i-22),:)=z(i,:);
        end
        %figure('Name', 'Upshifted Distribution of Single-Rates');
        %imagesc(temp);
        z = temp;
        upshift=1;
        izmax = izmax + 6;
      elseif izmax>21
        temp = zeros(28,8);
        for i=1:6
          temp((i+22),:)=z(i,:);
        end
        for i=7:1:28
          temp((i-6),:)=z(i,:);
        end
        %figure('Name', 'Downshifted Distribution of Single-Rates');
        %imagesc(temp);
        z = temp;
        downshift=1;
        izmax = izmax - 6;
      end
      z = z-zmin;
      zx = z(izmax(1),:)';
      zy = z(:,jzmax(1));
      x = (1:8)';
      y = (1:28)';
      fitx = fit(x,zx,'gauss1');
      fity = fit(y,zy,'gauss1','Exclude', zy<1000);

      %Evaluate fit
      stepx = 1/32.2857;
      stepy = 1/64.5714;
      xeval = (0:stepx:9.02);
      yeval = (0:stepy:28.01);
      zxeval= feval(fitx,xeval);
      zyeval= feval(fity,yeval);
    
      %Account for shift of x-component of the z-Matrix
      %  388 corresponds to  6  -> yeval( 388)= 5.9934
      %  389 corresponds to  7  -> yeval( 389)= 6.0089
      % 1421 corresponds to 22  -> yeval(1421)=21.9912
      % 1422 corresponds to 23  -> yeval(1422)=22.0066
      if downshift
        tempzyeval = zeros(length(zyeval),1);
        for i = 1:1421
          tempzyeval(i+388) = zyeval(i);
        end
        for i = 1422:1:length(zyeval)
          tempzyeval(i-1421) = zyeval(i);
        end
        zyeval = tempzyeval;
        %fprintf('%u: transaxial values are shifted by 6 detectors\r', k);
      elseif upshift
        tempzyeval = zeros(length(zyeval),1);
        for i = 1:388          
          tempzyeval(i+1421) = zyeval(i);
        end
        for i = 389:1:length(zyeval)
          tempzyeval(i-388) = zyeval(i);
        end
        zyeval = tempzyeval;
        %fprintf('%u: transaxial values are shifted by 6 detectors\r', k);
      end
      %figure('Name','Evaluated xFit');
      %plot(xeval,zxeval);
      %figure('Name','Evaluated yFit');
      %plot(yeval,zyeval);
    
      fitx1D(k) = xeval(find(zxeval==max(zxeval)));
      fity1D(k) = yeval(find(zyeval==max(zyeval)));
    end
  end
  % End Gaussian-Fit____________________________________
  %figure('Name','Evaluated Transaxial Fit');
  %plot(fity1D);
  
  %Add First and Last point do fity1D
  %timeList(timesteps+1) = tstop;
  %timeList(2:end+1)     = timeList;
  %timeList(1)           = 0;
  %fity1D(timesteps+1)   = 14.5;
  %fity1D(2:end+1)       = fity1D;
  %fity1D(1)             = 14.5;
  
  % ToDo:
  % Calculate and add First and Last Point to fitx1D!
  % Caluculation can be done using minlist

  % Flatten the fity to prepare for global fitting
  rounds = 1;
  flaty1D = fity1D;
  if backward
    for k = length(flaty1D):-1:2
      if flaty1D(k-1)<flaty1D(k)
        flaty1D(k-1) = flaty1D(k-1) + rounds*28;
      end
      if flaty1D(k-1)<flaty1D(k)
        rounds = rounds + 1;
        flaty1D(k-1) = flaty1D(k-1) + 28;
      end
    end
  else
    for k = 2:1:length(flaty1D)
      if flaty1D(k)<flaty1D(k-1)
        flaty1D(k) = flaty1D(k) + rounds*28;
      end
      if flaty1D(k)<flaty1D(k-1)
        rounds = rounds + 1;
        flaty1D(k) = flaty1D(k) + 28;
      end
    end
  end
  figure('Name','Evaluated Transaxial Fit flattened');
  plot(timeList,flaty1D);
  ylabel('Accumulated transaxial Bucket Number');
  xlabel('Time [ms]');

  % Global Fitting
  fitGlobalx = fit(timeList',fitx1D','poly3');
  fitGlobaly = fit(timeList',flaty1D','poly3');
  figure('Name','Global Fit X');
  plot(fitGlobalx,timeList,fitx1D');
  figure('Name','Global Fit Y');
  plot(fitGlobaly,timeList,flaty1D');
  %Evaluating global Fit
  fitGlobalxeval= feval(fitGlobalx,timeList);
  fitGlobalyeval= feval(fitGlobaly,timeList);
  %fitminyeval = feval(fitGlobaly,minlist);
  %figure();
  %plot(minlist,fitminyeval);
  %taking into account gravitation
  sinfity = flaty1D-fitGlobalyeval';
  plot(timeList,sinfity);
  xlabel('Time [ms]');
  fitGlobalSiny = fit(timeList',sinfity','sin2');
  figure('Name','Global Sin Fit Y');
  plot(fitGlobalSiny,timeList,sinfity');
  fitsinyeval= feval(fitGlobalSiny,timeList);
  figure();
  plot(sinfity-fitsinyeval');
  finalyfitflat = fitsinyeval'+fitGlobalyeval';
  figure();
  plot(timeList,finalyfitflat);
  % Unflatten fity
  rounds = 1;
  finalyfit = finalyfitflat;
  if backward
    for k = (length(finalyfit)-1):-1:1
      if finalyfit(k)>28
        finalyfit(k) = finalyfit(k) - rounds*28;
      end
      if finalyfit(k)>28
        rounds = rounds + 1;
        finalyfit(k) = finalyfit(k) - 28;
      end
    end
  else
    for k = 2:length(finalyfit)
      if finalyfit(k)>28
        finalyfit(k) = finalyfit(k) - rounds*28;
      end
      if finalyfit(k)>28
        rounds = rounds + 1;
        finalyfit(k) = finalyfit(k) - 28;
      end
    end    
  end
  figure();
  plot(timeList,finalyfit);

  figure;
  BlockSinglesMovie = zeros(28,8);
  F(timesteps) = struct('cdata',[],'colormap',[]);
  %G(time) = struct('cdata',[],'colormap',[]);
  for k=1:timesteps
    readPosition = 1;
    for i=1:8
      for j=1:28
        BlockSinglesMovie(j,i)=BlockSingles(k,readPosition);
        readPosition = readPosition + 1;
      end
    end
    contourf(BlockSinglesMovie);
    %surf(BlockSinglesMovie);
    %zlim([0 colorMax]);
    %zlabel('Bucket Singles Rate');
    %hold on;
    %imagesc(BlockSinglesMovie);
    %hold off;
    colormap(hot);
    legend(strcat('Time: ',num2str(int64(timeList(k))),' ms'),'Location','North');
    caxis([0 colorMax]);
    c = colorbar;
    c.Label.String = 'Bucket Singles Rate';
    xlabel('Ring Number');
    xlim([0 9]);
    ylabel('Single-Bucket Number');
    ylim([0 28]);
    %hold on;
    %plot(newfitx(k), fity1D(k), 'c*', 'MarkerSize', 20, 'LineWidth', 3)
    if k <= timesteps
      hold on;
      plot(fitGlobalxeval(k), finalyfit(k), 'b+', 'MarkerSize', 20, 'LineWidth', 3)
      hold off;
    end
    drawnow
    F(k) = getframe(gcf);
    %G(k) = getframe(gcf);
  end
    
  fig = figure('Name', 'Unrolled Rings of Single-Buckets');
  movie(fig,F,1,10)
  %movie(fig,G,1,24)
end

