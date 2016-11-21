function build_sino_static(filename, options)

Nbins=344;               % Number of radial bins
Nproj=252;               % Number of projections
Nplanes=4084;            % Number of 3D sinogram planes

paso_t=0.001;
sinoDim=Nbins*Nproj*Nplanes;

% Choose the file from an UI
%[filename, path]=uigetfile('*.*');
%cd(path);
% Get the Number of Counts (Ncs) from the size of the file
%cd('/home/andreas/data/PET_raw_data_20160603');
filesize=dir(filename);
Ncs=ceil(filesize.bytes/4);
fprintf('Estimated Number of Tags: \t %10.0f\r', Ncs);

fid=fopen(filename,'r');
dlist=fread(fid,[Ncs],'uint32');    
fclose(fid);
%cd('/home/andreas/code/PETreconstruction');

% -------------------------------------------------------
% Option 0: Output detailed information about all Tags
if options==0
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
        if     dlist(i) < 1073741824
            delays = delays + 1;
            %fprintf('Pos. %u\tat %s ms:\tDelay Event:\t%s\r',i,currentTime,num2str(dec2bin(dlist(i))));
        elseif dlist(i) < 2147483648
            prompts = prompts + 1;
            %fprintf('Pos. %u\tat %s ms:\tPrompt Event:\t%s\r',i,currentTime,num2str(dec2bin(dlist(i))));
        elseif dlist(i) < 2684354560
            timeTag = timeTag + 1;
            currentTime=num2str(dlist(i)-2^31);
            %fprintf('Pos. %u\tat %s ms:\tTime Tag:\t%s\r',i,currentTime,num2str(dec2bin(dlist(i))));
        elseif dlist(i) < 3221225472
            DtimeTag = DtimeTag + 1;
            %fprintf('Pos. %u\tat %s ms:\tDead-Time Tag:\t%s\r',i,currentTime,num2str(dec2bin(dlist(i))));
        elseif dlist(i) < 3758096384
            motionTag = motionTag + 1;
            %fprintf('Pos. %u\tat %s ms:\Motion Tag:\t%s\r',i,currentTime,num2str(dec2bin(dlist(i))));
        elseif dlist(i) < 4026531840
            patientTag = patientTag + 1;
            if     dlist(i)==3774877696
                fprintf('Pos. %u\tat %s ms:\tPatient Tag - Respiratory Trigger Gate on\r',i,currentTime);
            elseif dlist(i)==3774877697
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

% -------------------------------------------------------
% Option 1: Output prompts per second
elseif options==1
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
    
% -------------------------------------------------------
% Option 2: Output detailed dead-time information
elseif options==2
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
                fprintf('1. Lossy Node: %u \t lost events at %s [ms]\r', loss, Ttag);
            elseif typefield==6
                loss=bin2dec(binaryTag(13:32));
                LostEventsSecondLossyNode=LostEventsSecondLossyNode+loss;
                millionEvents=millionEvents+1;
                fprintf('2. Lossy Node: %u \t lost events at %s [ms]\r', loss, Ttag);
            %else
                %blocknum=bin2dec(binaryTag(4:13));
                %singles=bin2dec(binaryTag(14:32));
                %fprintf('Block: %u \tEvents: %u \t Time[ms]: %s\r', blocknum, singles, Ttag);
            end
        end
    end
    totalLoss=LostEventsSecondLossyNode+LostEventsFirstLossyNode;
    totalEvents=(millionEvents-1)*1000000;%1048575;
    fprintf('Total Dead-time marks: %u\r',D);
    fprintf('Number of lost events inserted by the "1. Lossy Node": %u\r', LostEventsFirstLossyNode);
    fprintf('Number of lost events inserted by the "2. Lossy Node": %u\r', LostEventsSecondLossyNode);
    fprintf('Total number of lost event packets: \t %u\r', totalLoss);
    fprintf('Total number of initial events: \t~%u\r', totalEvents);
    fprintf('Estimated Number of Tags: \t %u\r', Ncs);
    
% -------------------------------------------------------
% Option 3: Create Sinograms of the whole acquisition
elseif options==3
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
    
% -------------------------------------------------------
% Option 4: Cut the first X and the last Y ms of the acquisition
% and create Sinograms of the remaining events
else
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
    
    % Prompts between [001111...1] and [01111...1]
    ptag = cutdlist(find((cutdlist<2^31)&(cutdlist>=2^30)));
    % Randoms lower or equal [001111...1]
    rtag = cutdlist(find((cutdlist<2^30)));

    % TAGs (Time marks, dead-time tracking, fisio marks)
    Ttag = cutdlist(find((cutdlist>=2^31)&(cutdlist<2684354560)));
    Ttags=length(Ttag);
    time =Ttags*paso_t;
    clear Ttag;
    Dtag = cutdlist(find((cutdlist>=2684354560)&(cutdlist<3221225472)));
    Dtags=length(Dtag);
    clear Dtag;
    Ftag = cutdlist(find((cutdlist>=3758096384)&(cutdlist<3825205248)));
    Ftags=length(Ftag);
    clear Ftag;

    % remove preceding tag - necessary for ptag only, as the preceding
    % tag of rtag is already 0
    ptag = ptag-(2^30);
    % clear out-of-sinogram indices
    % -> bin-address of event packet has to be between 0 and sinoDim
    ptag = ptag((ptag<=sinoDim)&(ptag>0));
    rtag = rtag((rtag<=sinoDim)&(rtag>0));    

    % build sinograms
    sino=accumarray(ptag,1,[sinoDim,1])-accumarray(rtag,1,[sinoDim,1]);
    
    ptags = length(ptag);
    rtags = length(rtag);
    Ntags = ptags + rtags + Ttags + Dtags + Ftags;
    clear ptag;
    clear rtag;
    
    % write sinograms to file
    sinogramname = strcat('sinogram_static_', filename, '_cut.raw');
    fid=fopen(sinogramname,'w');
    fwrite(fid,uint16(sino),'uint16');
    fclose(fid);

    % Output for user
    fprintf('Finished reading list file\r\n');
    fprintf('Prompts\t\t:\t%s\r',num2str(ptags));
    fprintf('Randoms\t\t:\t%s\r',num2str(rtags));
    fprintf('ACQ time [s]\t:\t%s\r',num2str(time));
    fprintf('Total Number of Tags \t:\t%s\r\n',num2str(Ntags));
    
    if(Ntags<Ncs);
        %fprintf('Total number of Tags is smaller than "Ncs"!\n');
        fprintf('All Tags of the file were considered!\n');
    end
end

clear dlist


end
