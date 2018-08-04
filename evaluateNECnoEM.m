function evaluateNECnoEM()

halflifeFDG     = 6586.2;   % halflife of FDG in [s]
numberFiles     = 40;
results         = zeros(40,8);
scatterfraction = 0.13;
randomF         = 20./588;
% Description of Matrix Layout
%-1    -2        -3       -4       -5       -6     -7           -8
%-Time -Activity -Prompts -Randoms -Scatter -Trues -NECdahlborn -NECbailey

%_____________________________________________________
% (1) Time (= Start-Time from DICOM Header + 30 s)
% Start-Time is the time in seconds since 13:00
results(:,1)=[3342,3464,3591,3708,3833,3951,4056,4170, ...
              4287,4406,4531,4651,4768,4893,5010,5129, ...
              5279,5394,5502,5630,5764,5868,6002,6127, ...
              6240,6371,6502,6611,6739,6872,6989,7124, ...
              7254,7363,7600,7731,7862,8020,8170,8480];

%_____________________________________________________
% Activity
%syringe1
initAsyringe(1) = 358;   % initial Activity of syringe1 in [MBq]
initTimeS(1)    = 2862;  % Time of initA in [s] since 13:00 PET/MR time
%syringe2
initAsyringe(2) = 254.5; % initial Activity of syringe2 in [MBq]
initTimeS(2)    = 2932;  % Time of initA in [s] since 13:00 PET/MR time
%syringe3
initAsyringe(3) = 166.7; % initial Activity of syringe3 in [MBq]
initTimeS(3)    = 3075;  % Time of initA in [s] since 13:00 PET/MR time
%syringe4 -> replaces syringe 2 starting with measurement 35
initAsyringe(4) =  49.4; % initial Activity of syringe3 in [MBq]
initTimeS(4)    = 7486;  % Time of initA in [s] since 13:00 PET/MR time
%decay law: A(t) = A0 * 2^(-t/halflife)

% (2) Activity -> Activity of the respective Syringe
k=1;
for i=1:11
    for j=1:3
        results(k,2)=initAsyringe(j)*2^(-((results(k,1)-initTimeS(j))/halflifeFDG));
        k=k+1;
    end
end
results(35,2)=initAsyringe(4)*2^(-((results(35,1)-initTimeS(4))/halflifeFDG));
results(36,2)=initAsyringe(3)*2^(-((results(36,1)-initTimeS(3))/halflifeFDG));
results(37,2)=initAsyringe(1)*2^(-((results(37,1)-initTimeS(1))/halflifeFDG));
results(38,2)=initAsyringe(4)*2^(-((results(38,1)-initTimeS(4))/halflifeFDG));
results(39,2)=initAsyringe(3)*2^(-((results(39,1)-initTimeS(3))/halflifeFDG));
results(40,2)=initAsyringe(4)*2^(-((results(40,1)-initTimeS(4))/halflifeFDG));

%_____________________________________________________
% (3) Prompts (from DICOM Header)
results(:,3)=[416246366,443950899,460631286,418287427,448867315, ...
              441523088,420478014,453098783,424947184,423224465, ...
              457597087,408078185,426252369,462315050,391565494, ...
              429105085,463712448,375067314,432377577,464829262, ...
              359945140,435257078,466063178,350063968,438965244, ...
              467472673,332517942,443880408,468974410,319043638, ...
              448588076,470334761,309756237,453341464,136231063, ...
              290326540,459546074,129062311,276935335,122884516];

%_____________________________________________________
% (4) Randoms (from DICOM Header)
results(:,4)=[257510355,234251179,191090184,255126483,230301193, ...
              179288132,252978005,227544484,169170277,250813954, ...
              224078108,158893232,248210311,220611533,149150518, ...
              245641024,216763880,139475091,243468699,213305744, ...
              130782318,240814740,209754522,123949511,237787926, ...
              206185275,115242066,234372493,202745115,107924759, ...
              230652244,199281907,101844918,227284681, 26104969, ...
               92746380,222375601, 23988201, 85852664, 22053191];

%_____________________________________________________
% (5) Scatter
for i=1:numberFiles
    results(:,5)=(results(:,3)-results(:,4))*scatterfraction;
end

%_____________________________________________________
% (6) Trues (Prompts-Randoms-Scatter)
results(:,6)=results(:,3)-results(:,4)-results(:,5);


% (7) NEC (-> Dahlbom2005)
for i=1:numberFiles
    results(i,7)=(results(i,8)*results(i,8)) ...
        /(results(i,8)+2*randomF*results(i,6)+results(i,7));
end
% (8) NEC (-> Bailey2003)
for i=1:numberFiles
    results(i,8)=(results(i,5)*(results(i,8)/(results(i,8)+results(i,7)))) ...
        *(results(i,5)*(results(i,8)/(results(i,8)+results(i,7)))) ...
        /(results(i,5)+2*randomF*results(i,6));
end


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

% show NEC curve comparison
figure();
h = stem(sortResults(:,2), [sortResults(:,5),sortResults(:,6),sortResults(:,7), ...
    sortResults(:,8),sortResults(:,9),sortResults(:,10),sortResults(:,11), ...
    sortResults(:,12),sortResults(:,13),sortResults(:,14)],'filled', 'LineStyle', 'none');
xlabel('Activity [MBq]');
ylabel('Counts');
h(1).Marker = '*';
h(2).Marker = '*';
h(3).Marker = '*';
h(4).Marker = '*';
h(5).Marker = 'o';
h(6).Marker = 'o';
h(7).Marker = 'o';
h(8).Marker = 'o';
h(9).Marker = 'o';
h(10).Marker = 'o';
legend('Prompts','Randoms','Scatter','Trues','NECDahlbom0034','NECBailey0034', ...
    'NECDahlbom002','NECBailey002','NECDahlbom1','NECBailey1', ...
    'Location','northeast');
xlim([120,550]);
ylim([0,500000000])


% show NEC curve (used for Publication)
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
xlim([120,220]);
ylim([50000000,250000000]);

end
