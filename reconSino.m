function reconSino(filename)
% Reconstruction of 'filename'

% Basic Parameters of Siemens Biograph mMR
Nbins   = 344;         % Number of radial bins (NRAD)
Nproj   = 252;         % Number of projections (NANG)
Nslices = 127;         % Number of slices (2*Nrings -1)

SSRBSino = zeros(Nbins,Nproj,Nslices,'double');

fid  = fopen(filename,'r');
for i=1:Nslices
  Sino2D = fread(fid,[Nbins,Nproj],'float32');
  SIN2D = double(Sino2D);
  
  %for v=1:Nbins
  %  for w=1:Nproj
  %    if isnan(Sino2D(v,w))
  %      SIN2D(v,w) = 0.0;
  %    end
  %  end
  %end

  %SIN2D = imresize(Sino2D,2,'nearest');
  %for v=1:Nbins
  %  for w=1:Nproj
  %    if Sino2D(v,w)==0
  %      SIN2D(v,w) = 0.0;
  %    else
  %      SIN2D(v,w) = 5*log(Sino2D(v,w));
  %    end
  %  end
  %end
  SSRBSino(:,:,i) = SSRBSino(:,:,i) + SIN2D;
  %if i<50
  %  SSRBSino(:,:,i) = SSRBSino(:,:,i) + SIN2D/(0.1020*i+1.9365);
  %elseif i<76
  %  SSRBSino(:,:,i) = SSRBSino(:,:,i) + SIN2D/7;
  %else
  %  SSRBSino(:,:,i) = SSRBSino(:,:,i) + SIN2D/(-0.1020*i+14.6531);
  %end
end
fclose(fid);

theta = 180./Nproj;
for i=1:Nslices
  recon(:,:,i)=iradon(SSRBSino(:,:,i),theta,Nbins);
end

name = strcat('recon_', filename);
fid = fopen(name, 'w');
fwrite(fid,recon,'float32');
fclose(fid);
end