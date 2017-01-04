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
xlabel('Time [s]');
ylabel('Prompts per second');
%findpeaks(temporalPromptDistribution,'Annotate','extents','WidthReference','halfheight');

% Cut the first X and the last Y ms of acquisition
dlist = cutlmdata(dlist);

% Output detailed dead-time information and
% Get TimeTag of SingleBuckets
singleBucketTimes = deadTimeInfo(dlist);

% Show Bucket-Single rates
singleTime = length(singleBucketTimes);
showBucketSingles(dlist,singleTime);

% Create Sinograms
a=makeSino(cutdlist);

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
% Output detailed dead-time information
function [singleBucket] = deadTimeInfo(dlist)
  D=0;
  bucketRounds = 0;
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
        %singles=bin2dec(binaryTag(14:32));
        %fprintf('Block: %u\tSingles: %u\t Time[ms]: %s\r',blocknum,singles,Ttag);
        if blocknum == 223
          bucketRounds = bucketRounds + 1;
          singleBucket(bucketRounds) = str2num(Ttag);
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
function cutdlist = cutlmdata(dlist)
  % Read values for X and Y
  promptx='Enter time in ms you want first cut:';
  prompty='Enter time in ms you want second cut:';
  x=input(promptx);
  y=input(prompty);
  
  % Convert time in [ms] to the format of a Ttag and search for
  % the corresponding position of the Ttag in dlist 
  posx = find(dlist==(x+2^31));
  posy = find(dlist==(y+2^31));
    
  % Create the cut data list
  cutdlist = dlist(posx:posy);
end

% -------------------------------------------------------
% Show single-bucket rates
function showBucketSingles(dlist,time)
  D=0;
  millionEvents=0;
  BlockSingles = zeros(time,224);
  k=1;
  for i=1:length(dlist)
    if (dlist(i)>=2684354560)&&(dlist(i)<3221225472)
      % Dead-Time Tag
      D=D+1;
      binaryTag=num2str(dec2bin(dlist(i)));
      typefield=bin2dec(binaryTag(4:6));
      if typefield==7
        millionEvents=millionEvents+1;
      elseif typefield==6
        millionEvents=millionEvents+1;
      else
        blocknum=bin2dec(binaryTag(4:13));
        singles=bin2dec(binaryTag(14:32));
        BlockSingles(k,blocknum+1) = singles;
        if blocknum == 223
          k = k+1;
        end
      end
    end
  end
  totalEvents=millionEvents*1048575;
  fprintf('Total Dead-time marks: %u\r',D);
  fprintf('Total number of initial events: \t~%u\r', totalEvents);
  
  figure();
  v = [0,3000,6000,9000,12000,14000,16000,17000,18000,19000,20000];
  contourf(BlockSingles, v);
  colormap(hot);
  caxis([0 20000]);
  c = colorbar;
  c.Label.String = 'Bucket Singles Rate';
  ylabel('Ring Number');
  xlabel('Time [s]');
  
%______________________________________________________
% Gaussian-Fit
backward = 0;
for k = 1:time
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
    %nofitx(k) = izmax(1);
    %nofity(k) = jzmax(1);
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
    %smooth09x = fit(x,zx,'smoothingspline','SmoothingParam',0.9);
    %smooth09y = fit(y,zy,'smoothingspline','SmoothingParam',0.9,'Exclude', zy<1000);
    %figure('Name','xFit');
    %plot(fitx,x',zx');

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
    %fprintf('%u: x0=%u, y0=%u\r',k,fitx1D(k),fity1D(k));
  end
end

figure;
plot(fity1D);

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
figure;
plot(flaty1D);

%Fitting the global Values
xfitx = (1:time)';
yfity = (1:time)';
fitGlobalx = fit(xfitx,fitx1D','poly3');
fitGlobaly = fit(yfity,flaty1D','poly3');
figure('Name','Global Fit X');
plot(fitGlobalx,xfitx',fitx1D');
figure('Name','Global Fit Y');
plot(fitGlobaly,yfity',flaty1D');
%Evaluating global Fit
tevalx = (1:time);
tevaly = (1:time);
fitGlobalxeval= feval(fitGlobalx,tevalx);
fitGlobalyeval= feval(fitGlobaly,tevaly);

rounds = 1;
for k = 1:length(fitGlobalyeval)
  if fitGlobalyeval(k)>28
    fitGlobalyeval(k) = fitGlobalyeval(k) - rounds*28;
  end
  if fitGlobalyeval(k)>28
    rounds = rounds + 1;
    fitGlobalyeval(k) = fitGlobalyeval(k) - 28;
  end
end
  %cftool
  %______________________________________________________

  figure;
  BlockSinglesMovie = zeros(28,8);
  F(time) = struct('cdata',[],'colormap',[]);
  %G(time) = struct('cdata',[],'colormap',[]);
  for k=1:time
    readPosition = 1;
    for i=1:8
      for j=1:28
        BlockSinglesMovie(j,i)=BlockSingles(k,readPosition);
        readPosition = readPosition + 1;
      end
    end
    %v = [0,1000,2000,3000,4000,5000,6000,7000,8000,9000, ...
    %    10000,11000,12000,13000,14000,15000,16000,17000, ...
    %    18000,19000,20000];
    contourf(BlockSinglesMovie);
    %surf(BlockSinglesMovie);
    %zlim([0 22000]);
    %zlabel('Bucket Singles Rate');
    %hold on;
    %imagesc(BlockSinglesMovie);
    %hold off;
    colormap(hot);
    legend(strcat('Time: ',num2str(k*2),' s'),'Location','North');
    caxis([0 22000]);
    c = colorbar;
    c.Label.String = 'Bucket Singles Rate';
    xlabel('Ring Number');
    xlim([0 9])
    ylabel('Single-Bucket Number');
    ylim([0 28])
    %hold on;
    %plot(newfitx(k), fity1D(k), 'c*', 'MarkerSize', 20, 'LineWidth', 3)
    if k <= time
      hold on;
      plot(fitGlobalxeval(k), fity1D(k), 'b+', 'MarkerSize', 20, 'LineWidth', 3)
      hold off;
    end
    drawnow
    F(k) = getframe(gcf);
    %G(k) = getframe(gcf);
  end
    
  fig = figure('Name', 'Unrolled Rings of Single-Buckets');
  movie(fig,F,1,1)
  %movie(fig,G,1,24)
end

