function[x] = peak_voxel( TIMF, func_atlas );


% Extract data from one .mat file in each volnnnn subdirectory of the input subject directory.

TIMF_total = zeros(91,109,91,1200); %TIMF-total = image(x-length, y-length, z-length, timeseries length)

  %func_atlas = load_nii('atlas_funct.nii.gz'); % input atlas AAL116
  %TIMF = load_nii('rfMRI_REST1_LR.nii.gz'); %input signal
  TIMF_total = TIMF.img;
  

%%finding peak voxel

R = zeros(1,116); %116 = number of regions in the atlas that is used
l= 0;
for i = 1:116,
    S = [];
    S = find(func_atlas.img==i);
    [x,y,z] = ind2sub([91,109,91],S);
    LR(i) = length(x);
    M((l+1):(l+length(x)),1) = x;
    M((l+1):(l+length(x)),2) = y;
    M((l+1):(l+length(x)),3) = z;
    M((l+1):(l+length(x)),4) = zeros(length(x),1)+i;
    l=l+LR(i);
end
message=[ 'End of finding atlas coordination' ];
disp(message)
tic

L_ROIs = sum(LR(:));
Tseries = zeros(L_ROIs,1205);

    parfor j = 1:L_ROIs,


            timeseries(j,:) =  TIMF_total(M(j,1),M(j,2),M(j,3),:);
        end

%% Final Matrix
Tseries(:,1:4) = M;
Tseries(:, 5:1204) = timeseries;

%% Identifying the peak voxel of each region
peakvoxel = [];

    for i = 1:L_ROIs,
        Tseries(i, 1205) = (rms(Tseries(i, 5:1204)))^2  ;
    end
   


    L = 1;
    LR = LR;
    for j = 1:116,
        MAX(j) = max(Tseries (L:(L+LR(j)-1),1205));
        L = L + LR(j)
        N = [];
        N = find (Tseries(: , 1205)== MAX(j));

        peakvoxel(j,:)= Tseries(N,1:1205); %peak voxel with its coordination
    end
  
message=[ 'End of finding peak voxels' ];
disp(message)

   
	x = peakvoxel(:, 5:1204);%peakvoxel of each region
        
    
message=[ 'End of peak_voxel function' ];
disp(message)
end
