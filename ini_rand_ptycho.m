%%%% ini_rand_ptycho(str1,str2,os_rate, num_mask,coer,SUBIM, type_im,l_patch,x_c_p, y_c_p,bd_con, mask_type, c_l, delta, beta_1,beta_2, rho)  
%%% rand ptychograph initialization. 

%%% subim = size of each small patch
%%% bd_con = 0 ; 0 boundary condition 
%%%        = 1 ; periodic boundary condtion
%%%        = 2 ; not consider boundary
%%% l_patch  size of patch l_patch - by - l_patch
%%% x_c_p  x coordinates of center of patch


%%% mask_type = 0   iid mask
%             = 1   correlated mask
%             = 2   fresnel mask
 

function [IM,fft2_IM,fft2_IM_pois,NSR,phase_arg,nor_ptycho,nor_phase]=ini_rand_ptycho(str1,str2,os_rate,...
    num_mask,coer,SUBIM, type_im,l_patch,x_c_p, y_c_p,bd_con, mask_type, c_l, beta_1,beta_2, rho)        

IM1=imread(str1);
if length(size(IM1))>2
    IM1=rgb2gray(IM1);   %%% convert the image to grayscale. 
end
IM1=IM1(1:SUBIM,1:SUBIM);     %%% consider a sub image problem in order to save time.
IM1=double(IM1);

IM2=imread(str2);
if length(size(IM2))>2
    IM2=rgb2gray(IM2);   %%% convert the image to grayscale. 
end
IM2=IM2(1:SUBIM,1:SUBIM);     %%% consider a sub image problem in order to save time.
IM2=double(IM2);

switch type_im
    case 0    % real image
        IM=IM1;
    case 1    % str1+i*str1
        IM=IM1+1i*IM1;
    case 2
        IM=IM1+1i*IM2;   % construct str1+1i*str2.  
    case 3
        IM=IM1.*exp(2i*pi*rand(size(IM1))); % str1 * rand_phase
end


IM=IM/coer;        %%% coer to regulate poisson SNR. coer =400; SNR=3-4; coer=1 == no noise
%IM=ones(size(IM))*100/coer;
%IM=IM.*exp(i*rand(IM)); 







if mask_type == 0


switch num_mask             % 0  0 mask / Id mask
                            % 1  1 mask case
                            % 2  2 mask case 
                            % 3  1 and 1/2 mask case
    case 0
        phase_arg = zeros(l_patch,l_patch);
        mask = ones(l_patch,l_patch);
    case 1
        phase_arg=rand(l_patch,l_patch);  % argument of phase
        mask = exp(2i*pi*phase_arg);
    case 2
        phase_arg=rand(l_patch,l_patch,2);
        mask = exp(2i*pi*phase_arg);
    case 3
        mask = ones(l_patch,l_patch,2);
        mask(:,:,1) = exp(2i*pi*rand(l_patch,l_patch));
    case 4
        mask=ones(l_patch,l_patch,3);
        mask(:,:,1)=exp(2i*pi*rand(l_patch,l_patch));
        mask(:,:,2)=exp(2i*pi*rand(l_patch,l_patch));
end

elseif mask_type == 1    %% correlated mask please refer to ini_rand_ptycho_salt
    switch num_mask             % 0  0 mask / Id mask
                            % 1  1 mask case
                            % 2  2 mask case 
                            % 3  1 and 1/2 mask case
    case 0
        phase_arg = zeros(l_patch,l_patch);
        mask = ones(l_patch,l_patch);
    case 1
        phase_arg=rand(l_patch,l_patch);  
        mask = exp(2i*pi*phase_arg);
    case 2
        phase_arg=rand(l_patch,l_patch,2);
        mask = exp(2i*pi*phase_arg);
    case 3
        mask = ones(l_patch,l_patch,2);
        mask(:,:,1) = exp(2i*pi*rand(l_patch,l_patch));
    case 4
        mask=ones(l_patch,l_patch,3);
        mask(:,:,1)=exp(2i*pi*rand(l_patch,l_patch));
        mask(:,:,2)=exp(2i*pi*rand(l_patch,l_patch));
    end

    
    
elseif mask_type == 2   %%% fresnel mask
    switch num_mask             % 0  0 mask / Id mask
                            % 1  1 mask case
                            % 2  2 mask case 
                            % 3  1 and 1/2 mask case
    case 0
        phase_arg = zeros(l_patch,l_patch);
        mask = ones(l_patch,l_patch);
    case 1
        phase_arg= 1/2 * rho/l_patch *( (((1:l_patch)-beta_1).^2)'*ones(1,l_patch)+...
            (ones(l_patch,1)*((1:l_patch)-beta_2).^2)   );  % argument of phase
        mask = exp(2i*pi*phase_arg);
    case 2
        phase_arg=rand(l_patch,l_patch,2);
        mask = exp(2i*pi*phase_arg);
    case 3
        mask = ones(l_patch,l_patch,2);
        mask(:,:,1) = exp(2i*pi*rand(l_patch,l_patch));
    case 4
        mask=ones(l_patch,l_patch,3);
        mask(:,:,1)=exp(2i*pi*rand(l_patch,l_patch));
        mask(:,:,2)=exp(2i*pi*rand(l_patch,l_patch));
    end
    
    
    
    
    
    
    
end



p_c=0;
%%% calculate nor_factors 
IM_=ones(size(IM));
[Na,Nb]=size(IM_);
p_fft_IM=ptycho_rand_fft(IM_,os_rate,num_mask,mask,l_patch,x_c_p,y_c_p,bd_con);
nor_ptycho=ptycho_rand_ifft(p_fft_IM,Na,Nb,os_rate,num_mask,mask,l_patch,x_c_p,y_c_p,p_c);


phase_=ones(l_patch,l_patch);
fft_phase=pr_phase_fft(phase_,os_rate,num_mask,IM,l_patch,x_c_p,y_c_p); 
nor_phase=pr_phase_ifft(fft_phase,os_rate,num_mask,IM,l_patch,x_c_p,y_c_p);


%%% end calculate nor_factors



fft2_IM=ptycho_rand_fft(IM,os_rate,num_mask,mask,l_patch,x_c_p,y_c_p,bd_con);
fft2_IM_pois1=abs(fft2_IM).^2;
fft2_IM_pois=poissrnd(fft2_IM_pois1,size(fft2_IM));                     % poisson noise on sector
NSR=norm(abs(fft2_IM)-sqrt(fft2_IM_pois),'fro')/norm(abs(fft2_IM),'fro'); 
% NSR of data. agree with Chen.














end