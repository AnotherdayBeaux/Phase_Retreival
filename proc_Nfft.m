%%% SEE the heading of proc_ %%%



%%% str:        location of the image.


function [IM,fft2_IM,fft2_IM_pois,SNR,mask]=proc_Nfft(str,os_rate,num_mask,add_phase,sector_alpha,sector_beta,sector_mode,coer,SUBIM)        


IM=imread(str);
if length(size(IM))>2
    IM=rgb2gray(IM);   %%% convert the image to grayscale. 
end
[sub_x,sub_y]=size(IM);

IM=IM(1:floor(sub_x/SUBIM),1:floor(sub_y/SUBIM));     %%% consider a sub image problem in order to save time.
IM=double(IM);
IM=IM/coer;        %%% coer to regulate poisson SNR. coer =400; SNR=3-4;


if add_phase == 1
    switch sector_mode  %%%  sector_mode 0 general sector constrains; %
                        %%%              1 real  positivity;          %
                        %%%              2 complex positivity; 
                        %%%              3 non constraint 
        case 0
            IM=IM.*exp(1i*pi*((sector_alpha+sector_beta)*rand(size(IM))-sector_alpha));
              
        case 2
            IM=IM.*exp(1i*pi/2*rand(size(IM)));
        case 3
            IM=IM.*exp(1i*2*pi*rand(size(IM)));
    end
end

switch num_mask             % 0  0 mask / Id mask
                            % 1  1 mask case
                            % 2  2 mask case 
                            % 3  1 and 1/2 mask case
    case 0
        mask = ones(size(IM));
    case 1
        mask = exp(2i*pi*rand(size(IM)));
    case 2
        mask = exp(2i*pi*rand([size(IM),2]));
    case 3
        mask = ones([size(IM),2]);
        mask(:,:,1) = exp(2i*pi*rand(size(IM)));
end

fft2_IM=Nos_fft(IM,os_rate,num_mask,mask);
  
fft2_IM_pois1=abs(fft2_IM).^2;
fft2_IM_pois=poissrnd(fft2_IM_pois1,size(fft2_IM));                     % poisson noise on sector
SNR=norm(abs(fft2_IM).^2,'fro')/norm(abs(fft2_IM).^2-fft2_IM_pois,'fro');
end