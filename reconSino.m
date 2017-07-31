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
  SIN2D  = double(Sino2D);
  SSRBSino(:,:,i) = SSRBSino(:,:,i) + SIN2D;
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