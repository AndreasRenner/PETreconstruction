function build_sino_static(filename, options)

NRAD=336;                % Number of radial bins
NANG=336;                % Number of projections
Nsinos=559;              % Number of 3D sinogram planes

paso_t      =   0.001;
sino_dims   =   [NRAD,NANG,Nsinos];

% Choose the file from an UI
%[filename, path]=uigetfile('*.*');
%cd(path);
% Get the Number of Counts (Ncs) from the size of the file
cd('/home/andreas/data/PET_raw_data_20160603');
filesize=dir(filename);
Ncs=ceil(filesize.bytes/4);
fprintf('Number of Counts: %s\r', num2str(Ncs));

fid=fopen(filename,'r');
dlist=fread(fid,[Ncs],'uint32');    
fclose(fid);
cd('/home/andreas/code/PETreconstruction');


if options==1
    % Output prompts per second
    ptag=[];
    Ttag=[];
    Ftag=[];
    p=1;
    T=1;
    F=1;
    tmplen=0;
    for i=1:length(dlist)
        if (dlist(i)<(2^31))&&(dlist(i)>=(2^30))
            ptag(p)=dlist(i);
            p=p+1;
        elseif (dlist(i)>=(2^31))&&(dlist(i)<2684354560)
            Ttag(T)=dlist(i);
            if mod(Ttag(T)-(2^31),1000)==0
                fprintf('Prompts in s %s: %u\r', num2str((T-1)/1000), (length(ptag)-tmplen));
                tmplen=length(ptag);
            end
            T=T+1;
        end
    end
    clear Ftag;
elseif options==2
    % Output detailed dead-time information
    Dtag = dlist(find((dlist>=2684354560)&(dlist<3221225472)));
    fprintf('Total Dead-time marks: %s\r',num2str(length(Dtag)));
    DeadTime = zeros(length(Dtag));
    D=1;
    LostEvents=0;
    for i=1:length(dlist)
        if (dlist(i)>=(2^31))&&(dlist(i)<2684354560)
            Ttag=dlist(i);
            Ttag=num2str(Ttag-2^31);
        elseif (dlist(i)>=2684354560)&&(dlist(i)<3221225472)
            DeadTime(D)=dlist(i);
            D=D+1;
            binaryTag=num2str(dec2bin(dlist(i)));
            blocknum=bin2dec(binaryTag(4:13));
            singles=bin2dec(binaryTag(14:32));
            if blocknum==896||blocknum==768
                LostEvents=LostEvents+singles;
            end
            fprintf('Block: %u \tEvents: %u \t Time[ms]: %s\r', blocknum, singles, Ttag);
        end
    end
    fprintf('Total number of lost events: %u\r', LostEvents);
else
    % Create Sinograms
    % Matrix ops-------------------------------------------------------
    % Prompts between [001111...1] and [01111...1]
    ptag = dlist(find((dlist<2^31)&(dlist>=2^30)));
    % Randoms lower or equal [001111...1]
    rtag = dlist(find((dlist<2^30)));

    % TAGs (Time marks, dead-time tracking, fisio marks)
    Ttag = dlist(find((dlist>=2^31)&(dlist<2684354560)));
    Dtag = dlist(find((dlist>=2684354560)&(dlist<3221225472))); %>=101 00...0 <110 00...0
    Ftag = dlist(find((dlist>=3758096384)&(dlist<3825205248))); %>=111 00...0 <1110 0100 00...0

    % remove preceding tag if necessary
    ptag = ptag-(2^30);
    %Ttag = Ttag-(2^31);   % returns time in ms
    %find(~mod(Ttag, 1000))     % returns time in s
    % clear out-of-sinograms indices
    ptag = ptag((ptag<=prod(sino_dims))&(ptag>0));
    rtag = rtag((rtag<=prod(sino_dims))&(rtag>0));    

    % build sinograms
    sino=accumarray(ptag,1,[NRAD*NANG*Nsinos,1])-accumarray(rtag,1,[NRAD*NANG*Nsinos,1]);

    fid=fopen('sinogram_static.raw','w');
    fwrite(fid,uint16(sino),'uint16');
    fclose(fid);

    % Output for user
    fprintf('Finished reading list file\r\n');
    fprintf('Prompts\t\t:\t%s\r',num2str(length(ptag)));
    fprintf('Randoms\t\t:\t%s\r',num2str(length(rtag)));   
    fprintf('\r');
    fprintf('Time marks\t:\t%s\r',num2str(length(Ttag)));
    fprintf('Dead-time marks\t:\t%s\r',num2str(length(Dtag)));    
    %fprintf('Fisio marks\t:\t%s\r\n',num2str(length(Ftag)));   
    time=length(Ttag)*paso_t;
    Nevnts=length(ptag)+length(rtag)+length(Ttag)+length(Dtag)+length(Ftag);
    fprintf('ACQ time [s]\t:\t%s\r',num2str(time));
    fprintf('Total events\t:\t%s\r\n',num2str(Nevnts));
    if(Nevnts<Ncs);
        fprintf('Numer of "Total events" is smaller than "Ncs"!\n');
    end
end

clear dlist


end
