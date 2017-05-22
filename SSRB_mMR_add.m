function SSRB_mMR_add(filename, Nscans)

% Basic Parameters of Siemens Biograph mMR
Nbins   = 344;         % Number of radial bins (NRAD)
Nproj   = 252;         % Number of projections (NANG)
Nslices = 127;         % Number of slices (2*Nrings -1)

% SSRB Algorithm: creates a rebinned direct sinogram  
SSRBSino = zeros(Nbins,Nproj,Nslices,'double');

for i=1:Nscans
  name = strcat('sino_SSRB_',filename,num2str(i),'.raw');
  fid  = fopen(name,'r');
  for j=1:Nslices
    Sino2D = fread(fid,[Nbins,Nproj],'float32');
    SIN2D  = double(Sino2D);
    SSRBSino(:,:,j) = SSRBSino(:,:,j) + SIN2D;
  end
  fclose(fid);
end

name = strcat('sino_SSRB_', filename, '.raw');
fid = fopen(name,'w');
fwrite(fid,SSRBSino,'float32');
fclose(fid);

end

