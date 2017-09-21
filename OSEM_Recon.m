function readOut = OSEM_Recon(sino,filename)

% PROGRAM FOR OSEM RECONSTRUCTION
% GFN, UCM, JUN 2012. (Contact info: jacobo@nuclear.fis.ucm.es)
%
% PARAMETERS TO RECONSTRUCTION (MODIFY TO ADAPT TO YOUR CASE):
% SINOGRAM PARAMETERS
% Nslices = 64;        % Sinogram - Number of rings
Nslices  = 127;        % Direct planes - (2*NUMBER OF RINGS -1)
NRAD     = 344;        % Sinogram - Number of radial bins
NANG     = 252;        % Sinogram - Number of angles

% OSEM RECONSTRUCTION PARAMETERS
% Options: fact. of NANG=1 2 3 4 6 7 9 12 14 18 21 24 28 36 42 63 84 126 252
niter    = 24;         % NUMBER OF ITERATIONS
nsubsets = 12;         % NUMBER OF SUBSETS

% RELATED VARIABLES (DO NOT MODIFY) ------------------------
delta_ang = 180./NANG;

% PROGRAM FOR OSEM RECONSTRUCTION OF DIRECT SINOGRAMS (ART AND EMML
% ITERATIVE METHODS)

theta = 0:delta_ang:179.9;
EMML  = 1.0*ones(NRAD,NRAD,Nslices);   % INITIAL IMAGE (EM-ML)

for k=1:Nslices
  SINUNO = radon(ones(NRAD,NRAD),theta);    % CALIB FOR ART
  [rx, ry] = size(SINUNO);
  CAL = ones(rx,ry);                      % CALIB FOR EM-ML %WAS ones(ry,ry)
  DENOM_EMML = iradon(CAL,theta,'linear','None',1.0,NRAD);  % EM-ML
end

SN = zeros(rx,NANG,Nslices,'double');
m  = floor((rx-NRAD)/2-0.5);
for k=1:Nslices
  for j=1:NANG
    for i=1:NRAD
      SN(i+m,j,k) = sino(i,j,k);
      if (SN(i+m,j,k)<0);
        SN(i+m,j,k) = 0.0;
      end
    end
  end
end

for iter=1:niter
  disp(['Iteration ' num2str(iter) ' of ' num2str(niter)])
  gsize1=3;
  gsize2=3;
  
  for subiter=1:nsubsets
    vtheta2 = subiter:nsubsets:NANG;
    theta2  = (vtheta2-1)*delta_ang;
    fff     = fspecial3('gaussian',[gsize1 gsize1 gsize2]);
    EMML    = imfilter(EMML,fff,'replicate');
    disp(['Subiteration ' num2str(subiter)])
    for k=1:Nslices
      PROY_EMML(:,:,k) = radon(EMML(:,:,k),theta2); PROY_EMML(PROY_EMML==0)=1.0;
      CORR_EMML(:,:,k) = SN(:,vtheta2,k)./PROY_EMML(:,:,k);    %EM-ML
      NUM_EMML(:,:,k)  = iradon(CORR_EMML(:,:,k),theta2,'linear','None',1.0,NRAD);
      EMML(:,:,k)      = EMML(:,:,k).*(NUM_EMML(:,:,k)./DENOM_EMML);    % EM-ML
    end
  end
  
  % Output each iteration
  %name = strcat('osem_',num2str(nsubsets),'sub_',filename,'_iter-',num2str(iter),'.raw');
  %fid  = fopen(name,'w');
  %fwrite(fid,EMML,'float32');
  %fclose(fid);

  if iter==niter
    readOut = EMML;
    name = strcat('osem_',num2str(nsubsets),'sub_',filename,'_iter-',num2str(iter),'.raw');
    fid  = fopen(name,'w');
    fwrite(fid,EMML,'float32');
    fclose(fid);
  end
end

end