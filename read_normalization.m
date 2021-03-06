function read_normalization(filename, system)

%Normalization can be read for two different systems:
% - 1: Biograph64_TruePoint-1094
% - 2: Biograph_mMR

if system==1
    fprintf('Selected system is Biograph64_TruePoint\r');
%NORMALIZATION COMPONENTS DESCRIPTION:=
%number of normalization components:=7
%normalization component[1]:= geometric effects
%normalization component[2]:= crystal interference
%normalization component[3]:= crystal efficiencies
%normalization component[4]:= axial effects
%normalization component[5]:= paralyzing ring Dead-Time parameters
%normalization component[6]:= non-paralyzing ring Dead-Time parameters
%normalization component[7]:= TX crystal Dead-Time parameter
%size[1]:= {336,109}          {sino. proj. bins,sinogram planes}
%size[2]:= { 14,336}          {crystal number1 ,sino. proj. bins}
%size[3]:= {672, 55}          {crystal number2 ,ring number}
%size[4]:= {559}              {plane number}
%size[5]:= { 55}              {ring number}
%size[6]:= { 55}              {ring number}
%size[7]:= { 14}              {crystal number1}

    % System specific variables
    Nbins     = 336;
    Nslices   = 109;
    Ncrystal1 = 14;
    Ncrystal2 = 672;
    Nring     = 55;
    Nplanes   = 559;
    
    length    = Nbins*Nslices + Ncrystal1*Nbins + Ncrystal2*Nring + ...
                Nplanes + 2*Nring + Ncrystal1;
            
    % reading of all 7 components into a list
    %cd('/home/andreas/data/PET_raw_data_20160603');
    fid = fopen(filename,'r');
    normlist = fread(fid,[length],'*float32');
    fclose(fid);
    %cd('/home/andreas/code/PETreconstruction');
    
    % (1) Geometric effects
    norm1 = zeros(Nbins,Nslices);
    readPosition = 1;
    for j=1:Nslices
        for i=1:Nbins
            norm1(i,j) = normlist(readPosition);
            readPosition=readPosition+1;
        end
    end
    disp(readPosition)
    
    % (2) Crystal interference
    norm2 = zeros(Ncrystal1,Nbins);
    for j=1:Nbins
        for i=1:Ncrystal1
            norm2(i,j) = normlist(readPosition);
            readPosition=readPosition+1;
        end
    end
    disp(readPosition)
    
    % (3) Crystal efficiencies
    norm3 = zeros(Ncrystal2,Nring);
    for j=1:Nring
        for i=1:Ncrystal2
            norm3(i,j) = normlist(readPosition);
            readPosition=readPosition+1;
        end
    end
    disp(readPosition)
    
    % (4) Axial effects
    norm4 = zeros(Nplanes,1);
    for i=1:Nplanes
        norm4(i) = normlist(readPosition);
        readPosition=readPosition+1;
    end
    disp(readPosition);
    
    % (5) Paralyzing ring Dead-Time parameters
    norm5 = zeros(Nring,1);
    for i=1:Nring
        norm5(i) = normlist(readPosition);
        readPosition=readPosition+1;
    end
    disp(readPosition);
    
    % (6) Non-paralyzing ring Dead-Time parameters
    norm6 = zeros(Nring,1);
    for i=1:Nring
        norm6(i) = normlist(readPosition);
        readPosition=readPosition+1;
    end
    %disp(norm6);
    
    % (7) TX crystal Dead-Time parameter
    norm7 = zeros(Ncrystal1,1);
    for i=1:Ncrystal1
        norm7(i) = normlist(readPosition);
        readPosition=readPosition+1;
    end
    %disp(norm7);
    
elseif system==2
    fprintf('Selected system is Biograph_mMR\r');
%NORMALIZATION COMPONENTS DESCRIPTION:=
%number of normalization components:=8
%normalization component[1]:= geometric effects
%normalization component[2]:= crystal interference
%normalization component[3]:= crystal efficiencies
%normalization component[4]:= axial effects
%normalization component[5]:= paralyzing ring Dead-Time parameters
%normalization component[6]:= non-paralyzing ring Dead-Time parameters
%normalization component[7]:= TX crystal Dead-Time parameter
%normalization component[8]:= additional axial effect
%size[1]:= {344,127}          {sino. proj. bins,sinogram planes}
%size[2]:= {  9,344}          {crystal number1 ,sino. proj. bins}
%size[3]:= {504, 64}          {crystal number2 ,ring number}
%size[4]:= {837}              {plane number}
%size[5]:= { 64}              {ring number}
%size[6]:= { 64}              {ring number}
%size[7]:= {  9}              {crystal number1}
%size[7]:= {837}              {plane number}

    % System specific variables
    Nbins     = 344;
    Nslices   = 127;
    Ncrystal1 = 9;
    Ncrystal2 = 504;
    Nring     = 64;
    Nplanes   = 837;
    
    length    = Nbins*Nslices + Ncrystal1*Nbins + Ncrystal2*Nring + ...
                Nplanes + 2*Nring + Ncrystal1 + Nplanes;
            
    % reading of all 8 components into a list
    %cd('/home/andreas/data/PET_raw_data_20160603');
    fid = fopen(filename,'r');
    normlist = fread(fid,[length],'*float32');
    fclose(fid);
    %cd('/home/andreas/code/PETreconstruction');
    
    % (1) Geometric effects
    norm1 = zeros(Nbins,Nslices);
    readPosition = 1;
    for j=1:Nslices
        for i=1:Nbins
            norm1(i,j) = normlist(readPosition);
            readPosition=readPosition+1;
        end
    end
    disp(readPosition)
    
    % (2) Crystal interference
    norm2 = zeros(Ncrystal1,Nbins);
    for j=1:Nbins
        for i=1:Ncrystal1
            norm2(i,j) = normlist(readPosition);
            readPosition=readPosition+1;
        end
    end
    disp(readPosition)
    
    % (3) Crystal efficiencies
    norm3 = zeros(Ncrystal2,Nring);
    for j=1:Nring
        for i=1:Ncrystal2
            norm3(i,j) = normlist(readPosition);
            readPosition=readPosition+1;
        end
    end
    disp(readPosition)
    
    % (4) Axial effects
    norm4 = zeros(Nplanes,1);
    for i=1:Nplanes
        norm4(i) = normlist(readPosition);
        readPosition=readPosition+1;
    end
    disp(readPosition);
    
    % (5) Paralyzing ring Dead-Time parameters
    norm5 = zeros(Nring,1);
    for i=1:Nring
        norm5(i) = normlist(readPosition);
        readPosition=readPosition+1;
    end
    disp(readPosition);
    
    % (6) Non-paralyzing ring Dead-Time parameters
    norm6 = zeros(Nring,1);
    for i=1:Nring
        norm6(i) = normlist(readPosition);
        readPosition=readPosition+1;
    end
    %disp(norm6);
    
    % (7) TX crystal Dead-Time parameter
    norm7 = zeros(Ncrystal1,1);
    for i=1:Ncrystal1
        norm7(i) = normlist(readPosition);
        readPosition=readPosition+1;
    end
    %disp(norm7);
    
    % (8) Additional axial effect
    norm8 = zeros(Nplanes,1);
    for i=1:Nplanes
        norm8(i) = normlist(readPosition);
        readPosition=readPosition+1;
    end
    %disp(norm8);

else
    fprintf('Type "1" for Biograph_TruePoint or "2" for Biograph_mMR!\n');
   
end
