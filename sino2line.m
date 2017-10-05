function [line1,line2] = sino2line(sinoTotal,Nparts)

line1 = zeros(Nparts,252);
line2 = zeros(Nparts,252);

% -> create disk-shaped structuring element SE to preserve
%    round edges of BW object
SE = strel('disk',3);

% Define x-vector for evaluation of fit
xeval = (1:252);

faktor = 2/(344*252);

for i=1:Nparts
  
  %tic
  
  sino(:,:,:) = sinoTotal(i,:,:,:);
  a = sum(sum(sino(:,:,7)));
  sino(:,:,7) = imgaussfilt(sino(:,:,7),2);
  
  % calculate threshold for segmentation
  threshold = a*faktor;

  % convert this slice to binary image
  % maybe for high countrates the threshold has to be set to 1
  BW = im2bw(sino(:,:,7),threshold);
  %clear sino
  %figure()
  %imshow(BW)
  
  % Another option would be to use higher threshold and use
  % bwconvhull() to get boundaries
  % -> this would require a transformation around the middle line
  %    to achieve convex behaviour for eache side
  % -> this option has a problem if the curvature on the concarve
  %    side is stronger than on the convex side
  % -> maybe also have a look at boundary()
  
  % perform a morphological close operation on the image
  BW = imclose(BW,SE);  
  %figure();
  %imshow(BW)
  
  % reduce image using SE to get rid of rough edges
  BW = imerode(BW,SE);

  % select larges area of image
  BW = bwpropfilt(BW,'FilledArea',1,'largest');
  %figure();
  %imshow(BW)

  % restore initial size of image
  BW = imdilate(BW,SE);
  
  % extract edges from binary image
  BWedge = edge(BW,'Sobel');
  %figure();
  %imshow(BWedge)

  % extract and return 2 longest edges
  BW1 = bwpropfilt(BWedge,'Extent',1,'largest');
  BW2 = bwpropfilt(BWedge,'Extent',2,'largest');
  BW2 = logical(BW2-BW1);
  %figure();
  %imshow(BW1)
  %figure();
  %imshow(BW2);
  
  [y1,x1] = find(BW1);
  [y2,x2] = find(BW2);
  
  fit1 = fit(x1,y1,'poly4');
  fit2 = fit(x2,y2,'poly4');
  
  eval1 = feval(fit1,xeval);
  eval2 = feval(fit2,xeval);
    
  % calculate derivative do get split point
  deriv1 = differentiate(fit1,xeval);
  deriv2 = differentiate(fit2,xeval);
  splitF1 = find(abs(deriv1)==min(abs(deriv1)));
  splitF2 = find(abs(deriv2)==min(abs(deriv2)));
  
  % Rearange date to get sinogram-line of start and end of sinogram window
  % -> the upper line is representing tstart -> line 1
  % -> the lower line is representing tstop  -> line 2
  if (splitF1>202 && splitF2>202)
    % check which line is the lower one
    if eval1(splitF1)>eval2(splitF2)
      line1(i,:) = eval1(:)';
      line2(i,:) = eval2(:)';
    else
      % change number of eval (1->2;2->1)
      line1(i,:) = eval2(:)';
      line2(i,:) = eval1(:)';
    end
  elseif (splitF1<50 && splitF2<50)
    % check which line is the lower one
    if eval1(splitF1)<eval2(splitF2)
      line1(i,:) = eval1(:)';
      line2(i,:) = eval2(:)';
    else
      % change number of eval (1->2;2->1)
      line1(i,:) = eval2(:)';
      line2(i,:) = eval1(:)';
    end
  elseif (splitF1<50 && splitF2>202)
    % check which line is the lower one
    if eval1(splitF1)<eval2(splitF1)
      line1(i,:) = eval1(:)';
      line2(i,:) = eval2(:)';
    else
      % change number of eval (1->2;2->1)
      line1(i,:) = eval2(:)';
      line2(i,:) = eval1(:)';
    end
  elseif (splitF1>202 && splitF2<50)
    % check which line is the lower one
    if eval1(splitF1)>eval2(splitF1)
      line1(i,:) = eval1(:)';
      line2(i,:) = eval2(:)';
    else
      % change number of eval (1->2;2->1)
      line1(i,:) = eval2(:)';
      line2(i,:) = eval1(:)';
    end
    
  else
    % get split point in not fitted data
    split1 = find(x1==splitF1,1,'first');
    split2 = find(x2==splitF2,1,'first');
    % calculate distance from middle line to check wich part
    % of the overlap to discard (on which line)
    if abs(eval1(splitF1)-172)>abs(eval2(splitF2)-172)
      % make a first estimation of the extent of the sinogram-window
      %NbinPos = 172 + (eval2(splitF2)-172)/2;
      %if split2>126
      %  dif=find(abs(eval2-NbinPos)<1,1,'first')-find(abs(eval1-NbinPos)<1,1,'first'); 
      %else
      %  dif=find(abs(eval1-NbinPos)<1,1,'last')-find(abs(eval2-NbinPos)<1,1,'last');          
      %end
      %fprintf('Sinogram window is %u pixel\r',dif);
      % use the estimation to discard overlap-data from outter curve
      %temp  = ceil(dif/2);
      size1 = size(y1,1);
      size2 = size(y2,1);
      
      % get the lower line
      newy1(1:(split1-20)) = y1(1:(split1-20));
      newx1(1:(split1-20)) = x1(1:(split1-20));
      tempsize1 = size(newy1,2);
      for k=split2:1:size2
        newy1(tempsize1) = y2(k);
        newx1(tempsize1) = x2(k);
        tempsize1 = tempsize1 + 1;
      end
      
      % get the upper line
      newy2(1:split2) = y2(1:split2);
      newx2(1:split2) = x2(1:split2);
      tempsize2 = size(newy2,2);
      for k=(split1+20):size1
        newy2(tempsize2) = y1(k);
        newx2(tempsize2) = x1(k);
        tempsize2 = tempsize2 + 1;
      end
      
    else
      % make a first estimation of the extent of the sinogram-window
      %NbinPos = 172 + (eval1(splitF1)-172)/2;
      %if split1>126
      %  dif=find(abs(eval1-NbinPos)<1,1,'first')-find(abs(eval2-NbinPos)<1,1,'first'); 
      %else
      %  dif=find(abs(eval2-NbinPos)<1,1,'last')-find(abs(eval1-NbinPos)<1,1,'last');          
      %end        
      %fprintf('Sinogram window is %u pixel\r',dif);
      % use the estimation to discard overlap-data from outter curve
      %temp  = ceil(dif/2);
      size1 = size(y1,1);
      size2 = size(y2,1);
      
      % get the lower line
      newy1(1:(split2-20)) = y2(1:(split2-20));
      newx1(1:(split2-20)) = x2(1:(split2-20));
      tempsize1 = size(newy1,2);
      for k=split1:size1
        newy1(tempsize1) = y1(k);
        newx1(tempsize1) = x1(k);
        tempsize1 = tempsize1 + 1;
      end
      
      % get the upper line
      newy2(1:split1) = y1(1:split1);
      newx2(1:split1) = x1(1:split1);
      tempsize2 = size(newy2,2);
      for k=(split2+20):size2
        newy2(tempsize2) = y2(k);
        newx2(tempsize2) = x2(k);
        tempsize2 = tempsize2 + 1;
      end
      
    end
    
    % perform fit
    newfit1 = fit(newx1',newy1','poly4');
    newfit2 = fit(newx2',newy2','poly4');
    % evaluate fit
    eval1 = feval(newfit1,xeval);
    eval2 = feval(newfit2,xeval);    
    % return result
    line1(i,:) = eval1(:)';
    line2(i,:) = eval2(:)';
      
  end
  
  %toc
        
  %figure
  %imshow(BW);
  %hold on
  %plot(xeval,eval1)
  %hold on
  %plot(xeval,eval2)
  %hold off
  
  clear newx1
  clear newx2
  clear newy1
  clear newy2

end

end