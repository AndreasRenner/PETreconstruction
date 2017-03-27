function get_Pellet_Position_mMR(dlist,minlist)
% This function needs a already cut dlist
% -> start of dlist = minlist(1)

% Output detailed dead-time information and
% Get TimeTag of SingleBuckets
[singleBucketTimes,singleBuckets] = deadTimeInfo(dlist);

% Transform PET-acquisition-time into Transmission-scan-time
singleBucketTimes = singleBucketTimes - minlist(1);
fprintf('First Minimum of Prompts/s was at %u [ms]\n',minlist(1));

% Show Bucket-Single rates
showBucketSingles(singleBuckets,singleBucketTimes,minlist);
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
        %fprintf('Pos: %u\tBlock: %u\tSingles: %u\t Time[ms]: %s\r',i,blocknum,singles,Ttag);
        if blocknum == 223
          singleBucket(bucketRounds) = str2num(Ttag);
          bucketRounds = bucketRounds + 1;
        end
      end
    end
  end
  totalLoss=LostEventsSecondLossyNode+LostEventsFirstLossyNode;
  totalEvents=millionEvents*1048575;
  %GIM_loss_fraction = (totalEvents-LostEventsFirstLossyNode)/totalEvents;
  %PDR_loss_fraction = (totalEvents-LostEventsSecondLossyNode)/totalEvents;
  % 1048575 corresponds to the 20 bit ,,lost event tally field''
  fprintf('Total Dead-time marks: %u\r',D);
  fprintf('Number of lost events inserted by the "1. Lossy Node": %u\r', LostEventsFirstLossyNode);
  fprintf('Number of lost events inserted by the "2. Lossy Node": %u\r', LostEventsSecondLossyNode);
  fprintf('Total number of lost event packets: \t %u\r', totalLoss);
  fprintf('Total number of initial events: \t~%u\r', totalEvents);
end

% -------------------------------------------------------
% Show single-bucket rates
function showBucketSingles(BlockSingles,timeList,minlist)
  % >>timeList<< contains the exat time in [ms] of the dead-time tag
  % which is used to get information about Single-Buckets
  % >>minlist<< contains the time in [ms] of each minimum
  % in the ploted prompts/s over time -> corresponds to round-times
  %   0 time is choosen to be the fitted first minimum
  timesteps = length(timeList);
  % Transform PET-acquisition-time into Transmission-scan-time
  minlist = minlist-minlist(1);
  
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
      x = (1:8)';
      y = (1:28)';
      
      % Use row/line of zmax for 1-D fit
      %zx = z(izmax(1),:)';
      %zy = z(:,jzmax(1));
      %fitx = fit(x,zx,'gauss1');
      %fity = fit(y,zy,'gauss1','Exclude', zy<1000);
      
      % Use sum of rows/lines for 1-D fit
      rowsum  = sum(z,2);
      linesum = sum(z);
      fitx = fit(x,linesum','gauss1');
      fity = fit(y,rowsum,'gauss1','Exclude', rowsum<5000);

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

