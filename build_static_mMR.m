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
fprintf('Estimated Number of Tags: %s\r', num2str(Ncs));

fid=fopen(filename,'r');
dlist=fread(fid,[Ncs],'uint32');    
fclose(fid);
%cd('/home/andreas/code/PETreconstruction');
%LMwordcount=487595518;
%header=Ncs-LMwordcount;
%disp(header);

% -------------------------------------------------------
% Option 1: Output prompts per second
if options==1
    ptag=[];
    Ttag=[];
    p=1;
    T=1;
    tmplen=0;
    for i=1:length(dlist)
        if (dlist(i)<(2^31))&&(dlist(i)>=(2^30))
            ptag(p)=dlist(i);
            p=p+1;
        elseif (dlist(i)>=(2^31))&&(dlist(i)<2684354560)
            Ttag(T)=dlist(i);
            fprintf('Prompts in ms %u: %u\r', (T-1), (length(ptag)-tmplen));
            tmplen=length(ptag);
            T=T+1;
        end
    end
    
% -------------------------------------------------------
% Option 2: Output detailed dead-time information
elseif options==2
    Dtag=dlist(find((dlist>=2684354560)&(dlist<3221225472)));
    %DeadTime=zeros(length(Dtag));
    D=1;
    Ttag=[];
    LostEventsFirstLossyNode=0;
    LostEventsSecondLossyNode=0;
    millionEvents=0;
    for i=1:length(dlist)
        if (dlist(i)>=(2^31))&&(dlist(i)<2684354560)
            % Time Tag
            Ttag=dlist(i);
            Ttag=num2str(Ttag-2^31);
        elseif (dlist(i)>=2684354560)&&(dlist(i)<3221225472)
            % Dead-Time Tag
            %DeadTime(D)=dlist(i);
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
            else
                blocknum=bin2dec(binaryTag(4:13));
                singles=bin2dec(binaryTag(14:32));
                fprintf('Block: %u \tEvents: %u \t Time[ms]: %s\r', blocknum, singles, Ttag);
            end
        end
    end
    totalLoss=LostEventsSecondLossyNode+LostEventsFirstLossyNode;
    totalEvents=millionEvents*1048575;
    fprintf('Total Dead-time marks: %s\r',num2str(length(Dtag)));
    fprintf('Number of lost events inserted by the "1. Lossy Node": %u\r', LostEventsFirstLossyNode);
    fprintf('Number of lost events inserted by the "2. Lossy Node": %u\r', LostEventsSecondLossyNode);
    fprintf('Total number of lost event packets: \t %u\r', totalLoss);
    fprintf('Total number of initial events: \t %u\r', totalEvents);
    
% -------------------------------------------------------
% Option 3: Create Sinograms of the whole acquisition
elseif options==3
    % Prompts between [001111...1] and [01111...1]
    ptag = dlist(find((dlist<2^31)&(dlist>=2^30)));
    % Randoms lower or equal [001111...1]
    rtag = dlist(find((dlist<2^30)));

    % TAGs (Time marks, dead-time tracking, fisio marks)
    Ttag = dlist(find((dlist>=2^31)&(dlist<2684354560)));
    Dtag = dlist(find((dlist>=2684354560)&(dlist<3221225472))); %>=101 0..0 <110 0..0
    Ftag = dlist(find((dlist>=3758096384)&(dlist<3825205248))); %>=111 0..0 <1110 010..0

    % remove preceding tag if necessary
    ptag = ptag-(2^30);
    %Ttag = Ttag-(2^31);   % returns time in ms
    %find(~mod(Ttag, 1000))     % returns time in s
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
    time=length(Ttag)*paso_t;
    Ntags=length(ptag)+length(rtag)+length(Ttag)+length(Dtag)+length(Ftag);
    fprintf('ACQ time [s]\t:\t%s\r',num2str(time));
    fprintf('Total Number of Tags \t:\t%s\r\n',num2str(Ntags));
    if(Ntags<Ncs);
        %fprintf('Total number of Tags is smaller than "Ncs"!\n');
        fprintf('All Tags of the file were considered!\n');
    end
    
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
    Dtag = cutdlist(find((cutdlist>=2684354560)&(cutdlist<3221225472)));
    Ftag = cutdlist(find((cutdlist>=3758096384)&(cutdlist<3825205248)));

    % remove preceding tag - necessary for ptag only, as the preceding
    % tag of rtag is already 0
    ptag = ptag-(2^30);
    % clear out-of-sinogram indices
    % -> bin-address of event packet has to be between 0 and sinoDim
    ptag = ptag((ptag<=sinoDim)&(ptag>0));
    rtag = rtag((rtag<=sinoDim)&(rtag>0));    

    % build sinograms
    sino=accumarray(ptag,1,[sinoDim,1])-accumarray(rtag,1,[sinoDim,1]);

    % write sinograms to file
    sinogramname = strcat('sinogram_static_', filename, '_cut.raw');
    fid=fopen(sinogramname,'w');
    fwrite(fid,uint16(sino),'uint16');
    fclose(fid);

    % Output for user
    fprintf('Finished reading list file\r\n');
    fprintf('Prompts\t\t:\t%s\r',num2str(length(ptag)));
    fprintf('Randoms\t\t:\t%s\r',num2str(length(rtag)));  
    time=length(Ttag)*paso_t;
    Ntags=length(ptag)+length(rtag)+length(Ttag)+length(Dtag)+length(Ftag);
    fprintf('ACQ time [s]\t:\t%s\r',num2str(time));
    fprintf('Total Number of Tags \t:\t%s\r\n',num2str(Ntags));
    % Output lost events
    DeadTime=zeros(length(Dtag));
    D=1;
    LostEvents=0;
    for i=1:length(cutdlist)
        if (cutdlist(i)>=2684354560)&&(cutdlist(i)<3221225472)
            DeadTime(D)=cutdlist(i);
            D=D+1;
            binaryTag=num2str(dec2bin(cutdlist(i)));
            blocknum=bin2dec(binaryTag(4:13));
            singles=bin2dec(binaryTag(14:32));
            if blocknum==896||blocknum==768
                LostEvents=LostEvents+singles;
            end
        end
    end
    fprintf('Total number of lost events: %u\r\n', LostEvents);
    
    if(Ntags<Ncs);
        %fprintf('Total number of Tags is smaller than "Ncs"!\n');
        fprintf('All Tags of the file were considered!\n');
    end
end

clear dlist


end
