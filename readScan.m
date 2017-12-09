function [sinoTotal,maxIndex,tstart,tstop]=readScan(filename, ...
    Nscan,Nscans,Nparts,faktor,kumFaktor,direction, ...
    line1and2ref,maxIndexRef)

acqTime01Blank1 = 80170;
tpart = acqTime01Blank1/Nparts;

% Estimate Number of Tags from Filesize
file  = strcat(filename,'.bf');
size  = dir(file);
Ncs   = ceil(size.bytes/4);
fprintf('%s - Total Number of Tags: \t %10.0f\r',filename,Ncs);

Npos  = Ncs/Nscans;
Nread = ceil(Npos*faktor); % for 04HOT
fprintf('Number of Tags per read: %u\r\n',Nread);
fid   = fopen(file,'r');

offset = uint64(ceil(Npos*kumFaktor))*4;
fseek(fid,offset,'bof');

dlist = fread(fid,[Nread],'uint32');
  
% Get acquisition time in ms and output frequency of Tags
fprintf('\n');
fprintf('Information about Part %u:\r',Nscan);
acqTime = tagFrequency(dlist);
fprintf('Acquisition Time Part %u: %u [ms]\r',Nscan,acqTime);
  
name = strcat('minlistScan',num2str(Nscan),filename);
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

scanstart = minlist(1);
scanstop  = minlist(21);
% calculate conversionfactor for reference scan time
scantimef = (scanstop-scanstart)/acqTime01Blank1;

sinoTotal = zeros(Nparts,344,252,13);
maxIndex  = zeros(Nparts,1);
tstart    = zeros(Nparts,1);
tstop     = zeros(Nparts,1);
cut1 = scanstart;
tlen = (scanstop-scanstart)/Nparts;

for j=1:Nparts
  fprintf('Part %s:\r',num2str(j));
  %line1and2new = zeros(1,2,2523); %for Blank
  % tstart is only assigned a value if previous value is 0
  % otherwise it hase the value obtained from registration
  if ~tstart(j)
    tstart(j)= round(cut1);
    cut2 = cut1 + tlen;
  else
    cut2 = tstart(j) + tlen;
  end
  tstop(j) = round(cut2);
  dlistcut = cutlmdata(dlist,tstart(j),tstop(j));
  cut1 = cut2;
  sino = makeSinoRandomSubstraction(dlistcut);
  clear dlistcut;
  
  if ~line1and2ref
    [SSRBsino,maxIndex(j)] = SSRBmMR(sino,0);
    clear sino;
    sinoTotal(j,:,:,:) = SSRBsino;
  else
    %fprintf('Entering Sinogram-Window registration\r');
    if direction
      [SSRBsino,~] = SSRBmMR(sino,maxIndexRef(j));
      line1ref = line1and2ref(j,1,:);
      line2ref = line1and2ref(j,2,:);
    else
      [SSRBsino,~] = SSRBmMR(sino,maxIndexRef((Nparts+1)-j));
      line1ref = line1and2ref((Nparts+1)-j,1,:);
      line2ref = line1and2ref((Nparts+1)-j,2,:);
    end
    clear sino;
    % Check present sinogram window
    % calculate pixel2time faktor
    p2tfaktor = tpart/30; %line1ref(end);
    % get line1 and line2 of new scan
    %[line1and2new] = sino2line(SSRBsino,j,j); %forBlank
    [line1and2new] = sino2lineEM(SSRBsino,j,j);
    line1new(:,1) = line1and2new(:,1);
    line1new(:,2) = line1and2new(:,2);
    line2new(:,1) = line1and2new(:,1);
    line2new(:,2) = line1and2new(:,3);
    %difstart  = calculatedif(line1ref,line1and2new(1,1,:));
    %difstop   = calculatedif(line2ref,line1and2new(1,2,:));
    difstart  = calculatedifEM(line1ref,line1new);
    difstop   = calculatedifEM(line2ref,line2new);
    if abs(difstart-difstop)>20
      fprintf('Wrong line selected - I will try again\r');
      difstart  = calculatedifEM(line2ref,line1new);
      difstop   = calculatedifEM(line1ref,line2new);
      %difstart  = calculatedif(line2ref,line1and2new(1,1,:));
      %difstop   = calculatedif(line1ref,line1and2new(1,2,:));
    end
    clear line1new;
    clear line2new;
    if difstart<0
      if direction
        newstart = tstart(j) - abs(difstart)*p2tfaktor*scantimef;
      else
        newstart = tstart(j) + abs(difstart)*p2tfaktor*scantimef;
      end
    else
      if direction
        newstart = tstart(j) + difstart*p2tfaktor*scantimef;
      else
        newstart = tstart(j) - difstart*p2tfaktor*scantimef;
      end
    end
    if difstop <0
      if direction
        newstop = tstop(j) - abs(difstop)*p2tfaktor*scantimef;
      else
        newstop = tstop(j) + abs(difstop)*p2tfaktor*scantimef;
      end
    else
      if direction
        newstop = tstop(j) + difstop*p2tfaktor*scantimef;
      else
        newstop = tstop(j) - difstop*p2tfaktor*scantimef;
      end
    end
  
    if abs(difstart)<1 && abs(difstop)<1 % corresponds to 3% error; corresponds to 1 pixel
      if direction
        sinoTotal(j,:,:,:) = SSRBsino;
      else
        sinoTotal((Nparts+1)-j,:,:,:) = SSRBsino;
      end
    else
      clear SSRBsino;
      fprintf('Time start was: %7.0f\r',tstart(j));
      fprintf('Time stop  was: %7.0f\r',tstop(j));
      tstart(j)= newstart;
      tstop(j) = newstop;
      fprintf('Time start is:  %7.0f\r',tstart(j));
      fprintf('Time stop  is:  %7.0f\r',tstop(j));
      if j<Nparts
        tstart(j+1) = newstop;
      end
      dlistcut = cutlmdata(dlist,tstart(j),tstop(j));
      sino = makeSinoRandomSubstraction(dlistcut);
      clear dlistcut;
      
      if direction
        [SSRBsino,~] = SSRBmMR(sino,maxIndexRef(j));
        sinoTotal(j,:,:,:) = SSRBsino;
      else
        [SSRBsino,~] = SSRBmMR(sino,maxIndexRef((Nparts+1)-j));
        sinoTotal((Nparts+1)-j,:,:,:) = SSRBsino;
      end
      clear sino;
    end
  end
  fprintf('\r');
  
end

fclose(fid);

end

% *******************************************************
% ***         H E L P E R   F U N C T I O N S         ***
% *******************************************************

% -------------------------------------------------------
% calculate difference between present scan and reference scan
function dif = calculatedif(lineref,linenew)
splitF1 = lineref(end-1);
evalpnt = linenew(end);
%fprintf('Maximum Ref. is %f\r',splitF1);
%fprintf('Evaluation point is %f\r',evalpnt);

if splitF1>126
  pline1 = find(abs(lineref(1:end-2)-evalpnt)<0.11,1,'first')/10;
  pline2 = find(abs(linenew(1:end-2)-evalpnt)<0.11,1,'first')/10;
  if isempty(pline1)
    fprintf('Did not find anything for pline1\r');
    pline1 = find(abs(lineref(1:end-2)-evalpnt)<0.21,1,'first')/10;
  elseif isempty(pline2)
    fprintf('Did not find anything for pline2\r');
    pline2 = find(abs(linenew(1:end-2)-evalpnt)<0.21,1,'first')/10;
  end
else
  pline1 = find(abs(lineref(1:end-2)-evalpnt)<0.11,1,'last')/10;
  pline2 = find(abs(linenew(1:end-2)-evalpnt)<0.11,1,'last')/10;
  if isempty(pline1)
    fprintf('Did not find anything for pline1\r');
    pline1 = find(abs(lineref(1:end-2)-evalpnt)<0.21,1,'last')/10;
  elseif isempty(pline2)
    fprintf('Did not find anything for pline2\r');
    pline2 = find(abs(linenew(1:end-2)-evalpnt)<0.21,1,'last')/10;
  end
end


dif = pline1 - pline2;
fprintf('Difference is %f (%3.1f - %3.1f)\r',dif,pline1,pline2);
if dif>40
  fprintf('STOP\r');
end

end

% -------------------------------------------------------
% calculate difference between present scan and reference scan
function dif = calculatedifEM(lineref,linenew)
index = 1;
metricInit = 0;
for i=1:length(linenew(:,1))
  if lineref(linenew(i,1)*10+1)
    localdif = abs(linenew(i,2)-lineref(linenew(i,1)*10+1));
    if localdif>1
      metricInit = metricInit + sqrt(localdif);
    else
      metricInit = metricInit + localdif;
    end
    index = index + 1;
  end
end
metricInit = metricInit/index;

metricPlus = 0;
metricMinus= 0;
plus = 1;
minus= 1;

metricOld = metricInit;
while metricPlus<metricInit
  index = 1;
  for i=1:length(linenew(:,1))
    if ((i+plus)*10+1)<2522 && (i+plus)<length(linenew(:,1))
      localdif = abs(linenew(i,2)-lineref(linenew(i+plus,1)*10+1));
      if localdif>1
        metricPlus = metricPlus + sqrt(localdif);
      else
        metricPlus = metricPlus + localdif;
      end
      index = index + 1;
    end
  end
  metricPlus = metricPlus/index;
  if metricOld<metricPlus
    break
  end
  metricOld = metricPlus;
  plus = plus + 1;
end

if plus == 1
  metricOld = metricInit;
  while metricMinus<metricInit
    index = 1;
    for i=1:length(linenew(:,1))
      if (i-minus)>0
        localdif = abs(linenew(i,2)-lineref(linenew(i-minus,1)*10+1));
        if localdif>1
          metricMinus = metricMinus + sqrt(localdif);
        else
          metricMinus = metricMinus + localdif;
        end
        index = index + 1;
      end
    end
    metricMinus = metricMinus/index;
    if metricOld<metricMinus
      break
    end
    metricOld = metricMinus;
    minus = minus + 1;
  end
  dif = 0 - minus;
else
  dif = plus;
end

if plus==1 && minus ==1
  dif = 0;
end

fprintf('Difference is %f\r',dif);

end

% -------------------------------------------------------
% calculate difference between present scan and reference scan
function dif = calculatedifold(lineref,linenew)
splitF1 = lineref(2522);
splitF2 = linenew(2522);

% find optimal point to get difference
if splitF1>132 && splitF2>120
  point1 = (lineref(1)+linenew(1))/2;
  point2 = (lineref(splitF1*10)+linenew(splitF2*10))/2;
  difpnt = round((point1+point2)/2);
  pline1 = find(abs(lineref-difpnt)<0.11,1,'first')/10;
  pline2 = find(abs(linenew-difpnt)<0.11,1,'first')/10;
elseif splitF1>120 && splitF2>132
  point1 = (lineref(1)+linenew(1))/2;
  point2 = (lineref(splitF1*10)+linenew(splitF2*10))/2;
  difpnt = round((point1+point2)/2);
  pline1 = find(abs(lineref-difpnt)<0.11,1,'first')/10;
  pline2 = find(abs(linenew-difpnt)<0.11,1,'first')/10;
elseif splitF1>222 && splitF2<30
  point1 = (lineref(1)+linenew(1))/2;
  point2 = (lineref(2521)+linenew(2521))/2;
  difpnt = round((point1+point2)/2);
  pline1 = find(abs(lineref-difpnt)<0.11,1,'first')/10;
  pline2 = find(abs(linenew-difpnt)<0.11,1,'first')/10;
elseif splitF1<30 && splitF2>222
  point1 = (lineref(1)+linenew(1))/2;
  point2 = (lineref(2521)+linenew(2521))/2;
  difpnt = round((point1+point2)/2);
  pline1 = find(abs(lineref-difpnt)<0.11,1,'first')/10;
  pline2 = find(abs(linenew-difpnt)<0.11,1,'first')/10;   
else
  point1 = (lineref(2521)+linenew(2521))/2;
  point2 = (lineref(splitF1*10)+linenew(splitF2*10))/2;
  difpnt = round((point1+point2)/2);
  pline1 = find(abs(lineref-difpnt)<0.11,1,'last')/10;
  pline2 = find(abs(linenew-difpnt)<0.11,1,'last')/10;
end

if isempty(pline1)
  fprintf('Did not find anything for pline1\r');
  line(:) = lineref(1,1,1:2521);
  xeval2 = (0:0.1:252);
  figure();
  plot(xeval2,line);
elseif isempty(pline2)
  fprintf('Did not find anything for pline2\r');
  line(:) = linenew(1,1,1:2521);
  xeval2 = (0:0.1:252);
  figure();
  plot(xeval2,line);
end
dif = pline1 - pline2;
fprintf('Difference is %f (%3.1f - %3.1f)\r',dif,pline1,pline2);

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
if ~posx
  fprintf('Did not find tstart!\r');
end
posy = find(dlist==(int64(tstop)+2^31));
if ~posy
  fprintf('Did not find tstop!\r');
end
% Cut data list
dlist = dlist(posx:posy);  
end
