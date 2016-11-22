%Updated 23rd August by Amit Agrawal
%Null model for spatial indicators for the different transects
%For this moran's I code need to be in path https://in.mathworks.com/matlabcentral/fileexchange/13663-moran-s-i
clear all;
v = [1 2 3 4 5 6 7 8 9]; % to define the number of transect file number
for n=5 % for transect-5
    v1 = v(n)
    clearvars -except v n v1
    data=dlmread(strcat('ind_along83step_sliding_tran_',num2str(v1),'_reduced5.txt')); % Created thses files using the 'Real_Sp_Indicators.m' file
    Old_Spatial_Mean=data(:,2); % Colum second of real spatial indicator file contains the spatial mean (grass cover) 250x250 
    % defining empty vector for the all the values of spatial mean for a given transect
    spatial_mean = zeros(1,length(Old_Spatial_Mean)); 
    spatial_var = zeros(1,length(Old_Spatial_Mean)); 
    spatial_skew=zeros(1,length(Old_Spatial_Mean));
    spatial_corr=zeros(1,length(Old_Spatial_Mean));
    spatial_dft=zeros(1,length(Old_Spatial_Mean));
    % defining different possiblites of near neighbour for an element in the matrix
    W_raw=true(3); 
    W_raw(2,2)=0; W_raw(1,1)=0; W_raw(1,3)=0;
    W_raw(3,1)=0; W_raw(3,3)=0;
    norm=sum(W_raw(:)); 
    W=W_raw/norm;
    
% Simulate the 100 random matrices of a defined spatial mean (grass cover) but randomly filled
% vegetation with grass cover (1) and tree cover (0). The spatial mean was
% computed as average grass cover of 7.5 km x 7.5km sliding window over a
% transect.
% Following loop also calculates the mean and standard deviation of these
% 100 random matrix
% 100 randome matr
    for a=1:100; 
        for a0=1:length(Old_Spatial_Mean)
           f=floor(Old_Spatial_Mean(a0)*(250^2));
            idx = randperm(250*250);
            null_data = false(250,250); % create the matrix with full of zeros     
            null_data(idx(1:f)) = true; % fill the fraction of matrix randomly 
        
            null_width = length(null_data(1,:)); % width of null model
            null_length = length(null_data(:,1));% length of null model
    

                redu=5; % coarse grained hte matrix with 5 x 5
                x =zeros(1,(null_length-mod(null_length,redu))/redu); % defining the empty coarse grained matrix in x dimenesion
                y =zeros(1,(null_width-mod(null_width,redu))/redu);% defining the empty coarse grained matrix in y dimenesion
                local=zeros((null_length-mod(null_length,redu))/redu,(null_width-mod(null_width,redu))/redu); %definin the empty matrix for each spatial mean input
    
                for a1=1:redu:null_length-mod(null_length,redu) 
                    for b1=1:redu:null_width-mod(null_width,redu)
                        x=(a1+redu-1)/redu;
                        y=(b1+redu-1)/redu;
    
                        local(x,y) = sum(sum(null_data(a1:a1+redu-1,b1:b1+redu-1)))/redu^2; % coarse graining the matrix with redu = 5
            
                    end
                end
    
            
                spatial_mean(a0) = mean(local(:)); % compute the spatial mean of coarse grained matrix for a given spatial mean input from a real transect
                
                spatial_var(a0) = var(local(:));% compute the mean of spatial variance 
                
                spatial_skew(a0)=skewness(local(:));% compute the mean of spatial skewness
                
                mcor=moransI(null_data,W,'true'); % Keep this function in the path
                spatial_corr(a0)=mean(mcor(:));  % compute the mean of spatial correlation
                
                dft =abs(fft2(local)); %compute the root mean square of fast fourier transform of a matrix
                corner_local=length(local); % length of matrix
                range_opted=floor(length(local(1,:))*0.1); % matrix 10% of length chosen near to corner of matrix has been taken in account
                % defining the four corner submatrices in the matrix with
                % dimension is about 10% of the length of the matrix:
                % a_local, b_local, c_local and d_local
                before_corner=corner_local-range_opted+1; 
                a_local=dft(1:range_opted,1:range_opted);
                a_local(1)=[];
                b_local=dft( before_corner:corner_local,1:range_opted);
                b_local(length(b_local))=[];
                c_local=dft(before_corner:corner_local, before_corner:corner_local);
                c_local(length(c_local)^2)=[];
                d_local=dft(1:range_opted, before_corner:corner_local);
                d_local((length(d_local)-1)*length(d_local)+1)=[];

                abcd_local=[a_local;b_local;c_local;d_local]; % New matrix that is composed of four submatrices (dft)
                spatial_dft(a0) = mean(abcd_local(:)); % mean of the four submatrices that represent the spatial low frequency spectra
        end                

% form a matrix with that each row represent the indicaotrs value for the
% all the spatial mean from a given transect;  'a' is 100 times of simulation
        matrix_var(a,:)=spatial_var;
        matrix_mean(a,:)=spatial_mean;
        matrix_skew(a,:)=spatial_skew;
        matrix_corr(a,:)=spatial_corr;
        matrix_dft(a,:)=spatial_dft;   
    end

% compute the mean and standard deviation of indicators for the computes
% indicators from null matrix
    mean_matrix_mean=mean(matrix_mean);
    mean_matrix_var=mean(matrix_var);
    mean_matrix_skew=mean(matrix_skew);
    mean_matrix_corr=mean((matrix_corr));
    mean_matrix_dft=mean(matrix_dft);
    sd_mean=std(matrix_mean);
    sd_var=std(matrix_var);
    sd_skew=std(matrix_skew);
    sd_corr=std(matrix_corr);
    sd_dft=std(matrix_dft);
    
data=[mean_matrix_mean(:),sd_mean(:),mean_matrix_var(:),sd_var(:),mean_matrix_skew(:),sd_skew(:),mean_matrix_corr(:),sd_corr(:),mean_matrix_dft(:),sd_dft(:)]; % made a new matrix with mean inidcators value and corresponding standared deviation
dlmwrite(strcat('null_tran_83step_',num2str(v1),'_redu5.txt'),data,'delimiter','\t','-append');  % write the tab separated txt file for individual transects
end