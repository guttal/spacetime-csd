% Analysis of vegetation and MAP data
clear all;
tran_nu =[1 2 3 4 5 6 7 8 9]; % Transect Number identity vector
for x =5 % loop to write the file for indicators along rainfall gradient
    tran = tran_nu(x)
    trandata = dlmread(strcat('ind_along83step_sliding_tran_',num2str(tran),'_reduced5.txt'));
    nulldata = dlmread(strcat('null_tran_83step_',num2str(tran),'_redu5.txt'));
    map = xlsread(strcat('tran',num2str(tran),'_rain_2500m.xlsx'));
 
    scale_tran = trandata(:,1); % distance of sliding window over the transect
    veg_tran = trandata(:,2); % average grass cover of sliding window over the transect
    var_tran = trandata(:,3); % avgerage spatial variance of sliding window over the transect
    skew_tran = trandata(:,4); %average spatial skewness of sliding window over the transect
    cor_tran = trandata(:,5); %average correlation of sliding window over the transect
    dft_tran = trandata(:,6);%average spatial low frequency spectra of sliding window over the transect
    
    %average of different indicators computed for the 100
    %random matrix of given average grass cover along transect
    veg_null = nulldata(:,1);%average spatial grass cover of sliding window over the null transect
    var_null = nulldata(:,3);%average spatial variance of sliding window over the null transect
    skew_null = nulldata(:,5);%average spatial skewnesss of sliding window over the null transect
    cor_null = nulldata(:,7);%average spatial MoransI correlation of sliding window over the null transect
    dft_null = nulldata(:,9);%average spatial low frequency spectra of sliding window over the null transect
    %standard deviation for different indicators computed for the 100
    %random matrix of given average grass cover along transect
    sd_veg_null = nulldata(:,2);
    sd_var_null = nulldata(:,4);
    sd_skew_null = nulldata(:,6);
    sd_cor_null = nulldata(:,8);
    sd_dft_null = nulldata(:,10);
    
    redu = 3; % coarese grained the rainfall transect 3 pixexl x 3 pixel  = 7.5km x 7.5km
    g = map(:,1:3); %consdier 1-3 first coloums of rainfall transect    
    g_width=length(g(1,:));
    g_length=length(g(:,1));
    
   
    if (g_length>g_width) % calculate the number of sliding windonw with a step of 2.5km
        num_snaps = g_length-g_width+1 % number of sliding windows
        else num_snaps = g_width-g_length+1
    end
    
    for i1=1:num_snaps % for each sliding window calculate the avergate rainfall in the area of 7.5km x 7.5km
        if (g_length>g_width)  
            local = g(i1:i1+g_width-1,:);
            else local=g(:,i1:i1+g_length-1);
        end
          map_redu(i1) = mean(local(:));
    end
    
   
    [grad_map,idx] = sort(map_redu(:),1); % arrange the rainfall in ascending order = grad_map
    %here idx represent the series original index but now arranged for the ascending rainfall
    
    % using arranged idx for ascending rainfall in transects and null model
    % for computign grass cover and indicators along rainfall gradient
    veg_grad_map = veg_tran(idx); 
    var_grad_map = var_tran(idx);
    skew_grad_map = skew_tran(idx);
    corr_grad_map = cor_tran(idx);
    dft_grad_map = dft_tran(idx);
         
    nveg_grad_map = veg_null(idx);
    nvar_grad_map = var_null(idx);
    nskew_grad_map = skew_null(idx);
    ncorr_grad_map = cor_null(idx);
    ndft_grad_map = dft_null(idx);
    sd_nveg_grad_map = sd_veg_null(idx);
    sd_nvar_grad_map = sd_var_null(idx);
    sd_nskew_grad_map = sd_skew_null(idx);
    sd_ncorr_grad_map = sd_cor_null(idx);
    sd_ndft_grad_map = sd_dft_null(idx);
         
    data1 = [grad_map(12:47) veg_grad_map(12:47) var_grad_map(12:47) skew_grad_map(12:47) corr_grad_map(12:47) dft_grad_map(12:47)]; % consider 12 onwards
    datan = [nveg_grad_map(12:47) sd_nveg_grad_map(12:47) nvar_grad_map(12:47) sd_nvar_grad_map(12:47) nskew_grad_map(12:47) ...
    sd_nskew_grad_map(12:47) ncorr_grad_map(12:47) sd_ncorr_grad_map(12:47) ndft_grad_map(12:47) sd_ndft_grad_map(12:47)];


   dlmwrite(strcat('Indicator_tran',num2str(tran),'_rainfall.txt'),data1,'delimiter','\t','-append'); % write the transect files
   dlmwrite(strcat('null_indicators_tran',num2str(tran),'_rainfall.txt'),datan,'delimiter','\t','-append');% write the null transect file
end