% Updated on 23rd August 2016 by Amit Agrawal
% To compute the indicators over the transect with a slding window (7.5km x 7.5km) at 2.5km
% step
%%
clear all;
for transect_num=1:9
    transdata = xlsread(strcat('tran',num2str(transect_num),'veg_30m.xlsx'));
    
    
    tran_width = length(transdata(1,:)); % transect width
    tran_length = 10*floor(length(transdata(:,1))/10); %transect length
    
    if tran_length>tran_width   % Transesct 6,7,8 are not verically classified in the txt data
        tran_width=250; % set transect width = 7.5km
    else
        tran_length=250;
    end  
    %to assign the value zero or one for non-binary in matrices, non binary
    %is due to an absence of the valid data points
    for k=1:tran_length
        for l=1:tran_width
            if (transdata(k,l) < 0)
                if (rand<0.5)
                    transdata(k,l)=0;
                else
                    transdata(k,l)=1;
                end
            end
        end
    end
    
    redu_vector=[1 2 5 10]; % Vector of different coarse graining length
    
    clear data; 

    for r=1:4; %Run for loop for different coarse graining length
        redu=redu_vector(r);
    
        x = zeros(1,(tran_length-mod(tran_length,redu))/redu);
        y = zeros(1,(tran_width-mod(tran_width,redu))/redu);
        g = zeros((tran_length-mod(tran_length,redu))/redu,(tran_width-mod(tran_width,redu))/redu);

        for a=1:redu:tran_length-mod(tran_length,redu)
            for b=1:redu:tran_width-mod(tran_width,redu)
                x=(a+redu-1)/redu;
                y=(b+redu-1)/redu;
                g(x,y) = sum(sum(transdata(a:a+redu-1,b:b+redu-1)))/redu^2; % coarse grained matrix
            end
        end
        stp = [83 41 16 8] ; % step length required for the sliding window to move approximately 2.5 km
        step =stp(r); 
        g_width=length(g(1,:));
        g_length=length(g(:,1));
        if (g_length>g_width)
            num_snaps = length(1:step:g_length) % Number of sliding windows for a given coarse grained transect 
        else num_snaps = length(1:step:g_width)
        end

        W=true(3); 
        W(2,2)=0; W(1,1)=0; W(1,3)=0;
        W(3,1)=0; W(3,3)=0;
        norm=sum(W(:));
        W=W/norm;
        cor=zeros(1,num_snaps);
      
        s_windows = 1:step:step*num_snaps-g_width;

        for i1=1:length(s_windows)

            if (g_length>g_width)
                local = g(s_windows(i1):s_windows(i1)+g_width-1,:);
                
            else local=g(:,i1:s_windows(i1)+g_length-1);
            end
            
            imshow(local)
            spatial_scale(i1) = (i1-1)*30*redu*step; % Distance travelled by sliding window along a transect
            spatial_mean(i1) = mean(local(:));% average grass cover of  a sliding window
            spatial_var(i1) = var(local(:));% average spatial variance of a sliding window
            spatial_skew(i1) = skewness(local(:));% average sptial skewness of a sliding window
            mcor=moransI(local,W,'true'); % Morans I correlation
            spatial_corr(i1)=mean(mcor(:));% average spatial correlation of sliding window
            dft =abs(fft2(local)); % absoulute fast fourier transform of a sliding window

            corner_local=length(local); % length of a matrix
            range_opted=floor(length(local(1,:))*0.1);% matrix 10% of length chosen near to corner of matrix has been taken in account
            % defining the four corner submatrices in the matrix with
            % dimension is about 10% of the length of the matrix:
            % a_local, b_local, c_local and d_local
            before_corner=corner_local-range_opted+1;
            a_local=dft(1:range_opted,1:range_opted);
            a_local(1)=[];% to put the corner values void as the corner values so high peak that is not correspond to the dft?
            b_local=dft(before_corner:corner_local,1:range_opted);
            b_local(length(b_local))=[];
            c_local=dft(before_corner:corner_local, before_corner:corner_local);
            c_local(length(c_local)^2)=[];
            d_local=dft(1:range_opted, before_corner:corner_local);
            d_local((length(d_local)-1)*length(d_local)+1)=[];

           abcd_local=[a_local;b_local;c_local;d_local]; % New matrix that is composed of four submatrices (dft)
           
           spatial_dft(i1) = mean(abcd_local(:));% mean of the four submatrices that represent the spatial low frequency spectra
           
        end   
        data=[spatial_scale(:),spatial_mean(:),spatial_var(:),spatial_skew(:),spatial_corr(:),spatial_dft(:)];
        
        dlmwrite(strcat('ind_along83step_sliding_tran_',num2str(transect_num),'_reduced',num2str(redu),'.txt'),data,'delimiter','\t','-append');  % write the tab separated txt file for individual transects

        
    end
            
end