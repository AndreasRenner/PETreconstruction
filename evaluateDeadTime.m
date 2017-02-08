function evaluateDeadTime()

numberFiles     = 30;
results         = zeros(30,15);
scatterfraction = 0.13;
randomF         = 20./588;
% Description of Matrix Layout
%-1    -2       -3      -4       -5       -6       -7       -8
%-Time -Activity-A_Trans-A_Emiss -Prompts -Randoms -Scatter -Trues
%-9         -10         -11         -12        -13       -14
%-LostTally -TotalEvents-CorPrompts -CorDelays -CorTrues -LostEvents

%_____________________________________________________
% (1) Time (= Start-Time from DICOM Header + 30 s)
% Start-Time is the time in seconds since 14:00
results(:,1)=[1567,1709,1929,2123,2327,2561, ...
              2746,2905,3065,3242,3435,3601, ...
              3779,3954,4141,4299,4500,4683, ...
              4863,5043,5223,5393,5589,5760, ...
              5897,6119,6297,6420,6696,6874];

%_____________________________________________________
% Activity
halflifeFDG    = 6586.2;   % halflife of FDG in [s]
%bottle
initAbottle    = 41.3;     % initial Activity of bottle in [MBq]
initTimeB      = 789;      % Time of initA in [s] since 14:00 PET/MR time
%syringe1
initAsyringe(1)= 530.9;    % initial Activity of syringe1 in [MBq]
initTimeS(1)   = 1063;     % Time of initA in [s] since 14:00 PET/MR time
%syringe2
initAsyringe(2)= 380.2;    % initial Activity of syringe2 in [MBq]
initTimeS(2)   = 1196;     % Time of initA in [s] since 14:00 PET/MR time
%syringe3
initAsyringe(3)= 187.2;    % initial Activity of syringe3 in [MBq]
initTimeS(3)   = 1353;     % Time of initA in [s] since 14:00 PET/MR time
%decay law: A(t) = A0 * 2^(-t/halflife)

% (4) A_Emiss -> Activity of the bottle
for i=1:numberFiles
    results(i,4) = initAbottle*2^(-((results(i,1)-initTimeB)/halflifeFDG));
end

% (3) A_Trans -> Activity of the respective Syringe
k=1;
for i=1:10
    for j=1:3
        results(k,3)=initAsyringe(j)*2^(-((results(k,1)-initTimeS(j))/halflifeFDG));
        k=k+1;
    end
end

% (2) Activity-> Total activity
results(:,2)=results(:,3)+results(:,4);

%_____________________________________________________
% (5) Prompts (from DICOM Header)
results(:,5)=[396305312,412100834,458169481,398664857,415292681,464853078, ...
              401556827,418385629,466884956,403966560,421306680,469308732, ...
              406802872,425118272,471728658,409480360,429192051,455295696, ...
              412624637,433879039,427617234,415565085,438313778,404824607, ...
              418402903,442812119,377344940,420938831,447428825,352868298];

%_____________________________________________________
% (6) Randoms (from DICOM Header)
results(:,6)=[287150997,269762869,232609180,283969669,265965697,226957812, ...
              280280073,262233144,221862180,277248949,258784324,216556889, ...
              274037366,255222903,211450629,270470399,251560649,198009388, ...
              267214209,248084103,180242299,263761498,244525478,165010433, ...
              260293777,240906833,149211730,256735542,237219338,134660350];

%_____________________________________________________
% (7) Scatter
for i=1:numberFiles
    results(:,7)=(results(:,5)-results(:,6))*scatterfraction;
end

%_____________________________________________________
% (8) Trues (Prompts-Randoms-Scatter)
results(:,8)=results(:,5)-results(:,6)-results(:,7);

%_____________________________________________________
% (9,10) Dead-Time information (Lost Events and Total Events)
for filenum=1:numberFiles
    filename = strcat(int2str(filenum), '.IMA');
    filesize=dir(filename);
    Ncs=ceil(filesize.bytes/4);
    fprintf('Estimated Number of Tags in file %u: %u\r',filenum,Ncs);
    fid=fopen(filename,'r');
    dlist=fread(fid,[Ncs],'uint32');    
    fclose(fid);
       
    LostEventsFirstLossyNode=0;
    LostEventsSecondLossyNode=0;
    millionEvents=0;
    
    for i=1:length(dlist)
        if (dlist(i)>=2684354560)&&(dlist(i)<3221225472)
            % Dead-Time Tag
            binaryTag=num2str(dec2bin(dlist(i)));
            typefield=bin2dec(binaryTag(4:6));
            if typefield==7
                loss=bin2dec(binaryTag(13:32));
                LostEventsFirstLossyNode=LostEventsFirstLossyNode+loss;
                millionEvents=millionEvents+1;
            elseif typefield==6
                loss=bin2dec(binaryTag(13:32));
                LostEventsSecondLossyNode=LostEventsSecondLossyNode+loss;
                millionEvents=millionEvents+1;
            end
        end
    end
    
    clear dlist;
    
    results(filenum,9)=LostEventsSecondLossyNode+LostEventsFirstLossyNode;
    results(filenum,10)=millionEvents*1048575;
    
    % Output for user
    fprintf('Lost events "1. Lossy Node": %u\r',LostEventsFirstLossyNode);
    fprintf('Lost events "2. Lossy Node": %u\r',LostEventsSecondLossyNode);
    fprintf('Total lost event packets: \t %u\r',results(filenum,9));
    fprintf('Total initial events: \t %u\r\n',  results(filenum,10));

end

%_____________________________________________________
% Lost event tally correction [Jones2012]
% (11) Corrected Prompts
for i=1:numberFiles
    results(i,11)=results(i,5)+results(i,5)*results(i,9)/results(i,10);
end

% (12) Corrected Delays
for i=1:numberFiles
    results(i,12)=results(i,6)+results(i,6)*results(i,9)/results(i,10);
end

% (13) Corrected Trues
results(:,13)=results(:,11)-results(:,12);

% (14) Lost Events
results(:,14)=results(:,10)-results(:,5)-results(:,6);

% (15) NEC
for i=1:numberFiles
    results(i,15)=(results(i,8)*results(i,8)) ...
        /(results(i,8)+2*randomF*results(i,6)+results(i,7));
end

% (16) Relative amount of Corrected Randoms
%for i=1:numberFiles
%    results(i,16)=(results(i,5)+results(i,6))*0.15;
%end

fid=fopen('resultsnew.dat','w');
fwrite(fid,results,'real*4');
fid=fopen('resultsnew.dat','r');
results = fread(fid,[30,15],'real*4');
fclose(fid);

sortResults = sortrows(results,2);

figure();
h = stem(sortResults(:,2), [sortResults(:,5),sortResults(:,6), ...
    sortResults(:,7),sortResults(:,8),sortResults(:,9), ...
    (sortResults(:,9)-sortResults(:,8)-sortResults(:,6)-sortResults(:,5)), ...
    sortResults(:,10),sortResults(:,11),sortResults(:,12)], ...
    'filled', 'LineStyle', 'none');
h(1).Marker = 'o';
h(2).Marker = 's';
h(3).Marker = 'd';
h(4).Marker = '+';
h(5).Marker = 'o';
h(6).Marker = '*';
h(7).Marker = 'x';
h(8).Marker = 'x';
h(9).Marker = 'x';
h(4).Color = 'k';
h(5).Color = 'k';
h(6).Color = 'k';
h(7).Color = h(1).Color;
h(8).Color = h(2).Color;
h(9).Color = h(3).Color;
xlabel('Activity [MBq]');
ylabel('Counts');
legend('Prompts','Randoms','Trues','Lost Events','Total Events', ...
    'Other losses','Cor. Prompt','Cor. Randoms','Cor. Trues', ...
    'Location','northwest');
xlim([120,550]);

% show NEC curve
figure();
h = stem(sortResults(:,2), [sortResults(:,5),sortResults(:,6),sortResults(:,7), ...
    sortResults(:,8),sortResults(:,15)],'filled', 'LineStyle', 'none');
xlabel('Activity [MBq]');
ylabel('Counts');
h(1).Marker = '*';
h(2).Marker = '*';
h(3).Marker = '*';
h(4).Marker = '*';
h(5).Marker = 'o';
legend('Prompts','Randoms','Scatter','Trues','NEC', ...
    'Location','northeast');
xlim([120,550]);
ylim([0,500000000])

% show NEC curve
figure();
h = plot(sortResults(:,2), [sortResults(:,6),sortResults(:,8), ...
    sortResults(:,15)]);
xlabel('Activity [MBq]');
ylabel('Counts');
h(1).Marker = '*';
h(2).Marker = '+';
h(3).Marker = 'o';
h(1).LineStyle = '--';
h(2).LineStyle = ':';
legend('Randoms','Trues','NEC', ...
    'Location','northwest');
xlim([120,410]);
ylim([90000000,290000000]);

% show results vs. corrected results for prompt/random/true
%figure();
%k = stem(sortResults(:,2), [sortResults(:,5),sortResults(:,6),sortResults(:,7), ...
%    sortResults(:,10),sortResults(:,11),sortResults(:,12)], ...
%    'filled', 'LineStyle', 'none');
%k(1).Marker = 'square';
%k(2).Marker = 'square';
%k(3).Marker = 'square';
%k(4).Marker = '+';
%k(5).Marker = '+';
%k(6).Marker = '+';
%k(4).Color = k(1).Color;
%k(5).Color = k(2).Color;
%k(6).Color = k(3).Color;
%xlabel('Activity [MBq]');
%ylabel('Counts');
%legend('Prompts','Randoms','Trues','Prompts Corr','Randoms Corr','Trues Corr', ...
%   'Location','northwest');
%ylim([100000000,700000000])

%remove Syringe1 from data
%result12=results;
%for i=1:3:28
%    result12(i,:)=0;
%end

end
