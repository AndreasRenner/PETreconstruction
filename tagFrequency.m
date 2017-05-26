function duration = tagFrequency(dlist)
% Get frequency of all individual Tags

% Event Packets:  0XXX ...
prompts    = 0; % 01XX ...
delays     = 0; % 00XX ...
% Tag Packets:    1XXX ...
timeTag    = 0; % 100X ...
DtimeTag   = 0; % 101X ...
motionTag  = 0; % 1100 ...
patientTag = 0; % 1110 ...
controlTag = 0; % 1111 ...
%currentTime= 0;

fprintf('****************************************************************\r');
fprintf('START Detailed Information about Frequency of individual Tags:\r');
%fprintf('Special Control/Patient Tags:\r');
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
    %currentTime=num2str(dlist(i)-2^31);
    %fprintf('Pos %u\tat %s ms:\tTime Tag:\t%s\r',i,currentTime,num2str(dec2bin(dlist(i))));
  elseif   dlist(i) < 3221225472
    DtimeTag = DtimeTag + 1;
    %fprintf('Pos %u\tat %s ms:\tDead-Time Tag:\t%s\r',i,currentTime,num2str(dec2bin(dlist(i))));
  elseif   dlist(i) < 3758096384
    motionTag = motionTag + 1;
    %fprintf('Pos %u\tat %s ms:\Motion Tag:\t%s\r',i,currentTime,num2str(dec2bin(dlist(i))));
  elseif   dlist(i) < 4026531840
    patientTag = patientTag + 1;
    %if     dlist(i)== 3774877696
    %  specialTagCounter = specialTagCounter + 1;
    %  fprintf('Pos %u\tat %s ms:\tPatientTag - Respiratory Trigger Gate on\r',i,currentTime);
    %elseif dlist(i)== 3774877697
    %  specialTagCounter = specialTagCounter + 1;
    %  fprintf('Pos %u\tat %s ms:\tPatientTag - Respiratory Trigger Gate off\r',i,currentTime);
    %end
    %fprintf('Pos %u\tat %s ms:\tPatient Tag:\t%s\r',i,currentTime,num2str(dec2bin(dlist(i))));
  else
    controlTag = controlTag + 1;
    %if     dlist(i)==4294901760
    %  fprintf('Pos %u\tat %s ms:\tControlTag - Start of Acquisition\r',i,currentTime);
    %elseif dlist(i)==4286545920
    %  fprintf('Pos %u\tat %s ms:\tControlTag - Time Synchronization with MR\r',i,currentTime);
    %end
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
fprintf('****************************************************************\r');

% Return scan-duration in [ms]
duration = timeTag;
end
