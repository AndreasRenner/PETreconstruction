function reconSinoRatio(SSRB_Blanc,SSRB_Trans)
% Reconstruct ratio between Blanc- and Transmission-Scan

%load('SSRB_Blanc.mat');
%load('SSRB_Trans.mat');
SSRB_Blanc=smooth3(SSRB_Blanc,'gaussian',[7 7 5],2);
SSRB_Trans=smooth3(SSRB_Trans,'gaussian',[7 7 5],2);
%SSRB_Trans=smooth3(SSRB_Trans,'gaussian',[5 5 5],1);
%SSRB_Blanc=smooth3(SSRB_Blanc,'gaussian',[5 5 5],1);
%SSRB_Trans=smooth3(SSRB_Trans,'gaussian',[3 3 3],0.42466);
%SSRB_Blanc=smooth3(SSRB_Blanc,'gaussian',[3 3 3],0.42466);
for u=1:size(SSRB_Blanc,1)
  for v=1:size(SSRB_Blanc,2)
    for w=1:size(SSRB_Blanc,3)
      if SSRB_Trans(u,v,w)<0.1 % SSRB_Blanc(u,v,w)<4 &&
        SSRB_Ratio(u,v,w)=log(SSRB_Blanc(u,v,w));
      else
        SSRB_Ratio(u,v,w)=log(SSRB_Blanc(u,v,w)./SSRB_Trans(u,v,w));
      end
    end
  end
end

newName = 'sino_SSRB_ratio04RQ_TEST_Blankonly2.raw';
fid = fopen(newName,'w');
fwrite(fid,SSRB_Ratio,'float32');
fclose(fid);

theta = 180./size(SSRB_Blanc,2);
for i=1:size(SSRB_Blanc,3)
  recon(:,:,i)=iradon(SSRB_Ratio(:,:,i),theta);
end

fid = fopen('reconRatio04RQ_TEST.raw', 'w');
fwrite(fid,recon,'float32');
fclose(fid);
end