function build_static_mMR(filename)
% Basic Parameters of Siemens Biograph mMR
Nbins   = 344;         % Number of radial bins
Nproj   = 252;         % Number of projections
Nplanes = 4084;        % Number of 3D sinogram planes

sinoDim = Nbins * Nproj * Nplanes;

% Choose the file from an UI
%[filename, path]=uigetfile('*.*');
%cd(path);
% Get the Number of Counts (Ncs) from the size of the file
%cd('/home/andreas/data/PET_raw_data_20160603');

% Estimate Number of Tags from Filesize
filesize = dir(filename);
Ncs      = ceil(filesize.bytes/4);
clear filesize;
fprintf('Estimated Number of Tags: \t %10.0f\r', Ncs);

% Read File into a List
fid   = fopen(filename,'r');
dlist = fread(fid,[Ncs],'uint32');    
fclose(fid);
clear fid;
%cd('/home/andreas/code/PETreconstruction');

% Get acquisition time in ms and output frequency of Tags
acqTime = tagFrequency(dlist);
fprintf('Acquisition time: %u [ms]', acqTime);

% Output prompts per second
p_s(dlist);

% Output detailed dead-time information
deadtimeinfo(dlist);

% Cut the first X and the last Y ms of acquisition
cutdlist = cutlmdata(dlist);

% Create Sinograms
makesino(cutdlist);

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
  for i=1:length(dlist)
    if       dlist(i) < 1073741824
      delays = delays + 1;
      %fprintf('Pos. %u\tat %s ms:\tDelay Event:\t%s\r',i,currentTime,num2str(dec2bin(dlist(i))));
    elseif   dlist(i) < 2147483648
      prompts = prompts + 1;
      %fprintf('Pos. %u\tat %s ms:\tPrompt Event:\t%s\r',i,currentTime,num2str(dec2bin(dlist(i))));
    elseif   dlist(i) < 2684354560
      timeTag = timeTag + 1;
      currentTime=num2str(dlist(i)-2^31);
      %fprintf('Pos. %u\tat %s ms:\tTime Tag:\t%s\r',i,currentTime,num2str(dec2bin(dlist(i))));
    elseif   dlist(i) < 3221225472
      DtimeTag = DtimeTag + 1;
      %fprintf('Pos. %u\tat %s ms:\tDead-Time Tag:\t%s\r',i,currentTime,num2str(dec2bin(dlist(i))));
    elseif   dlist(i) < 3758096384
      motionTag = motionTag + 1;
      %fprintf('Pos. %u\tat %s ms:\Motion Tag:\t%s\r',i,currentTime,num2str(dec2bin(dlist(i))));
    elseif   dlist(i) < 4026531840
      patientTag = patientTag + 1;
      if     dlist(i)== 3774877696
        fprintf('Pos. %u\tat %s ms:\tPatient Tag - Respiratory Trigger Gate on\r',i,currentTime);
      elseif dlist(i)== 3774877697
        fprintf('Pos. %u\tat %s ms:\tPatient Tag - Respiratory Trigger Gate off\r',i,currentTime);
      end
      %fprintf('Pos. %u\tat %s ms:\tPatient Tag:\t%s\r',i,currentTime,num2str(dec2bin(dlist(i))));
    else
      controlTag = controlTag + 1;
      if     dlist(i)==4294901760
        fprintf('Pos. %u\tat %s ms:\tControl Tag - Start of Acquisition\r',i,currentTime);
      elseif dlist(i)==4286545920
        fprintf('Pos. %u\tat %s ms:\tControl Tag - Time Synchronization with MR\r',i,currentTime);
      end
      %fprintf('Pos. %u\tat %s ms:\tControl Tag:\t%s\r',i,currentTime,num2str(dec2bin(dlist(i))));     
    end
  end
  events = prompts + delays;
  tags   = timeTag + DtimeTag + motionTag + patientTag + controlTag;
  fprintf('Total number of Event Packets: \t %10.0f\r', events);
  fprintf('Total number of Tag Packets:   \t %10.0f\r\n', tags);
  fprintf('Prompt Events:  \t %10.0f\r', prompts);
  fprintf('Delay Events:   \t %10.0f\r', delays);
  fprintf('Time Tags:      \t %10.0f\r', timeTag);
  fprintf('Dead-Time Tags: \t %10.0f\r', DtimeTag);
  fprintf('Motion Tags:    \t %10.0f\r', motionTag);
  fprintf('Patient Tags:   \t %10.0f\r', patientTag);
  fprintf('Control Tags:   \t %10.0f\r', controlTag);
  
  duration = timeTag;
end

% -------------------------------------------------------
% Get prompts per second
function p_s(dlist)
  p=0;
  T=0;
  tmplen=0;
  for i=1:length(dlist)
    if (dlist(i)<(2^31))&&(dlist(i)>=(2^30))
      p=p+1;
    elseif (dlist(i)>=(2^31))&&(dlist(i)<2684354560)
      Ttag=dlist(i);
      if mod(Ttag-(2^31),1000)==0
        fprintf('Prompts in s %u: %u\r', T/1000, (p-tmplen));
        tmplen=p;
      end
      T=T+1;
    end
  end
end
    
% -------------------------------------------------------
% Output detailed dead-time information
function deadtimeinfo(dlist)
  D=0;
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
        fprintf('1. Lossy Node: %u\tlost events at %s [ms]\r',loss,Ttag);
      elseif typefield==6
        loss=bin2dec(binaryTag(13:32));
        LostEventsSecondLossyNode=LostEventsSecondLossyNode+loss;
        millionEvents=millionEvents+1;
        fprintf('2. Lossy Node: %u\tlost events at %s [ms]\r',loss,Ttag);
      else
        blocknum=bin2dec(binaryTag(4:13));
        singles=bin2dec(binaryTag(14:32));
        fprintf('Block: %u\tSingles: %u\t Time[ms]: %s\r',blocknum,singles,Ttag);
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
  fprintf('Estimated Number of Tags: \t %u\r', Ncs);
end

% -------------------------------------------------------
% Create Sinograms of the whole acquisition
function makesino(dlist)
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
% Cut the first X and the last Y ms of acquisition
function cutdlist = cutlmdata(dlist)
  % Read values for X and Y
  promptx='Enter time in ms you want to cut at the beginning:';
  prompty='Enter time in ms you want to cut at the end:';
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
% Option 5: Show block singles
elseif options==5
  %time = 499; %BlancScan
  time = 88; %TransmissionScan
  %time = 29; %Dead-Time Scans
  %timef= 482; %BlancScan
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
  clear dlist;
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
backward = 1;
for k = 1:time
  z = zeros(28,8);
  readPosition = 1;
  for i=1:8
    for j=1:28
      z(j,i) = BlockSingles(k,readPosition);
      readPosition = readPosition + 1;
    end
  end
  figure('Name', 'Initial Distribution of Single-Rates');
  imagesc(z);
  zmin = min(min(z));
  zmax = max(max(z));
  downshift=0;
  upshift=0;
  if (zmax-zmin)<1000
    fprintf('%u: Peak is too small for fitting!\r',k);
    tempstart = k;
  else
    tstart = tempstart+1;
    tstop  = k;
    [izmax,jzmax] = find(z==zmax);
    nofitx(k) = izmax(1);
    nofity(k) = jzmax(1);
    if k==tstart+1 && jzmax>6
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
fprintf('Start: %u, stop: %u \r',tstart,tstop);

figure;
plot(fity1D);
% Flatten the fity to prepare for global fitting
rounds = 1;
flaty1D = fity1D;
if backward
  for k = (length(flaty1D)-2):-1:2
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

% cut ends where transmission source is leaving helix
for i = tstart:1:(tstop-2)
  newfity(i-tstart+1) = flaty1D(i);
end
figure('Name', 'Cutted Fit Y');
plot(newfity);
for i = tstart:1:(tstop-2)
  newfitx(i-tstart+1)=fitx1D(i);
end
figure('Name', 'Cutted Fit X');
plot(newfitx);

%Fitting the global Values
xnewfitx = (tstart:(tstop-2))';
ynewfity = (tstart:(tstop-2))';
fitGlobalx = fit(xnewfitx,newfitx','poly3');
fitGlobaly = fit(ynewfity,newfity','poly3');
figure('Name','Global Fit X');
plot(fitGlobalx,xnewfitx',newfitx');
figure('Name','Global Fit Y');
plot(fitGlobaly,ynewfity',newfity');
%Evaluating global Fit
tevalx = (0:tstop);
tevaly = (0:tstop);
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
        BlockSinglesMovie(j,i) = BlockSingles(k,readPosition);
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
    if k <= tstop
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

clear dlist

end
