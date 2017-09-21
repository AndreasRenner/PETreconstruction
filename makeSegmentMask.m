function makeSegmentMask(x,y,segments)
% Input coordinates of central axis of the head coil in recon-space
% to create Segment-Masks

% Basic Parameters of Siemens Biograph mMR
Nbins   = 344;         % Number of radial bins (NRAD)
Nproj   = 344;         % Number of projections (NANG)
Nslices = 127;         % Number of slices (2*Nrings -1)

segMask1 = zeros(Nbins,Nproj,Nslices);
segMask2 = ones(Nbins,Nproj,Nslices);

for k=1:Nslices
  for i=x:Nbins
    for j=y:Nproj
      segMask1(i,j,k)=1;
      segMask2(i,j,k)=0;
    end
  end
  for i=1:(x-1)
    for j=1:(y-1)
      segMask1(i,j,k)=1;
      segMask2(i,j,k)=0;
    end
  end
end
if segments==4
elseif segments==8
  for k=1:Nslices
    for i=x:Nbins
      jmax = (y+i-x);
      if jmax>Nproj
        jmax = Nproj;
      end
      for j=y:jmax
        jmin = (j-2*(j-y)-1);
        if jmin<1
          jmin = 1;
        end
        segMask1(i,j,k)=0;
        segMask2(i,j,k)=1;
        segMask1(i,jmin,k)=1;
        segMask2(i,jmin,k)=0;
      end
    end
    for i=1:(x-1)
      jmax = (y+x-i-1);
      if jmax>Nproj
        jmax = Nproj;
      end
      for j=1:jmax
        jmin = (j-2*(j-y)-1);
        if jmin<1
          jmin = 1;
        elseif jmin>(y-1)
          jmin = (y-1);
        end
        segMask1(i,j,k)=1;
        segMask2(i,j,k)=0;
        segMask1(i,jmin,k)=0;
        segMask2(i,jmin,k)=1;
      end
    end
  end
else
  fprintf('Wrong Number of segments');
end

% Write segment mask to file
nameSegMask1 = 'Mask1.raw';
nameSegMask2 = 'Mask2.raw';
fid1 = fopen(nameSegMask1,'w');
fid2 = fopen(nameSegMask2,'w');
fwrite(fid1,segMask1,'float32');
fwrite(fid2,segMask2,'float32');
fclose(fid1);
fclose(fid2);

end
