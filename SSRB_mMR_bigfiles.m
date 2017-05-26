function SSRB_mMR_bigfiles(filename, fileending, Nscans)

% Estimate Number of Tags from Filesize
file  = strcat(filename, fileending);
size  = dir(file);
Ncs   = ceil(size.bytes/4);
fprintf('%s - Total Number of Tags: \t %10.0f\r',filename,Ncs);

Npos  = Ncs/Nscans;
Nread = Npos;
%Nread = ceil(Npos*1.1);
fprintf('Number of Tags per read: %u\r\n', Nread);
start = 0;
%start = 154938+118585+110162+155863+87207+124791; %97032+157230+96259
%tstart=[55380,177360,297590,422600,540230,658250,768720,879910,1006230,1115500];
%tstop=[140530,259620,382720,506210,625560,740220,848090,971760,1087890,1197180];
%faktor=[1.05,1.05,1.05,1,0.98,0.9,1.08,0.9,1,1];
fid   = fopen(file,'r');

%offset= uint64(ceil(Npos*1.1) + 3*ceil(Npos*1.05) + Npos + ceil(Npos*0.98))*4;
%status= fseek(fid,offset,'bof');
%pos   = ftell(fid);

for i=1:Nscans
  %pos   = ftell(fid);
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
  prompt = 'Enter start time in ms:';
  tstart = input(prompt);
  clear prompt;
  prompt = 'Enter stop time in ms:';
  tstop  = input(prompt);
  %clear prompt;

  %prompt = 'Enter Nread Faktor:';
  %faktor = input(prompt);
  Nread  = Npos;
  %Nread  = ceil(Npos*faktor(i));

  dlist  = cutlmdata(dlist,tstart,tstop);

  name   = strcat(filename, num2str(i));
  makeSino(dlist,name,0);
  clear dlist;
  SSRB_mMR(name);

  start  = start+acqTime;
end
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
