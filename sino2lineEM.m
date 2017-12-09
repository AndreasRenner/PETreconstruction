function [line1and2] = sino2lineEM(sinoTotal,Nparts,Nstart)

% -> create disk-shaped structuring element SE to preserve
%    round edges of BW object
SEs = strel('disk',2); % for blank
%SEs = strel('disk',4); % for emission
SEb = strel('disk',3);
%SEb = strel('disk',4);

%faktor = 1.9/(344*252); % for blank
faktor = 1.6/(344*252); % for emission

for i=Nstart:Nparts
  if Nparts==Nstart
    sino = sinoTotal;
    name = strcat('BWedge_Scan2_',num2str(i),'.raw');
    i=1;
  else
    sino(:,:,:) = sinoTotal(i,:,:,:);
    name = strcat('BWedge_Scan1_',num2str(i),'.raw');
  end
  a = sum(sum(sino(:,:,7)));
  
  sino = smooth3(sino,'gaussian',[5 5 5],1.5);  
  sino(:,:,7) = imgaussfilt(sino(:,:,7),2);
  
  % calculate threshold for segmentation
  threshold = a*faktor;

  % convert this slice to binary image
  if threshold<0.95
    BW = im2bw(sino(:,:,7),threshold);
  else
    BW = im2bw(sino(:,:,7),0.95);
  end
  %clear sino
  %figure()
  %imshow(BW)
  
  % perform a morphological close operation on the image
  BW = imclose(BW,SEb);  
  %figure();
  %imshow(BW)

  % reduce image using SE to get rid of rough edges
  BW = imerode(BW,SEb);
  %figure();
  %imshow(BW)
  
  % select areas larger than 50 pixel
  BW = bwareaopen(BW,50);
  
  % restore initial size of image
  BW = imdilate(BW,SEs);
  %figure();
  %imshow(BW)
  % perform a morphological close operation on the image
  
  BW = imclose(BW,SEb); %for emission
  
  %figure();
  %imshow(BW)

  index = 1;
  controlindex = 1;
  for j=1:252
    point1 = find(BW(:,j),1,'first');
    point2 = find(BW(:,j),1,'last');
    if point1
      y1(index) = point1;
      y2(index) = point2;
      x1(index) = j;
      index = index + 1;
    else
      fprintf('We have a hole in the edges detected (at j=%i).\r',j);
      controlindex = controlindex + 1;
    end
  end
  
  line1and2 = zeros(length(y1),3);
  line1and2(:,1) = x1(:);
  line1and2(:,2) = y1(:);
  line1and2(:,3) = y2(:);  

  fid = fopen(name,'w');
  fwrite(fid,BW,'float32');
  fclose(fid);
end

end