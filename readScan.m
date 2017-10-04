function [sinoTotal,maxIndex]=readScan(filename,Nscan,Nscans,Nparts,...
    faktor,offset)

% Estimate Number of Tags from Filesize
file  = strcat(filename,'.bf');
size  = dir(file);
Ncs   = ceil(size.bytes/4);
fprintf('%s - Total Number of Tags: \t %10.0f\r',filename,Ncs);

Npos  = Ncs/Nscans;
Nread = ceil(Npos*faktor); % for 04HOT
fprintf('Number of Tags per read: %u\r\n',Nread);
fid   = fopen(file,'r');

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

tstart = minlist(1);
tstop = minlist(21);

sinoTotal = zeros(Nparts,344,252,13);
maxIndex  = zeros(Nparts,1);
cut1 = tstart;
tlen = (tstop-tstart)/Nparts;
for j=1:Nparts
  cut2 = cut1 + tlen;
  dlistcut = cutlmdata(dlist,round(cut1),round(cut2));
  cut1 = cut2;
  sino = makeSinoRandomSubstraction(dlistcut);
  clear dlistcut;
  
  [SSRBsino,maxIndex(j)] = SSRBmMR(sino);
  clear sino;
  sinoTotal(j,:,:,:) = SSRBsino;
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
