function makeSino(dlist,filename)
% Create Sinograms

% Basic Parameters of Siemens Biograph mMR
Nbins   = 344;         % Number of radial bins (NRAD)
Nproj   = 252;         % Number of projections (NANG)
Nplanes = 4084;        % Number of 3D sinogram planes (Nsinos)
sinoDim = Nbins * Nproj * Nplanes; % 354 033 792

% Prompts between [001111...1] and [01111...1]
ptag = dlist((dlist<2^31)&(dlist>=2^30));
ptag = ptag-(2^30);
% clear out-of-sinograms indices
% sinoDim equal to 1 0101 0001 1010 0010 0000 1000 0000
ptag = ptag((ptag<=prod(sinoDim))&(ptag>0));
% build sinograms
sino = accumarray(ptag,1,[sinoDim,1]);
clear ptag;

% Randoms lower or equal [001111...1]
rtag = dlist((dlist<2^30));
% clear out-of-sinograms indices
rtag = rtag((rtag<=prod(sinoDim))&(rtag>0));
% substract sinograms
sino = sino-accumarray(rtag,1,[sinoDim,1]);
clear rtag;

% write sinograms to file
sinogramname = strcat('sino_', filename, '.raw');
fid=fopen(sinogramname,'w');
fwrite(fid,uint16(sino),'uint16');
fclose(fid);
end
