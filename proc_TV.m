%%% combination of 'add_noise'+'add_phase_mask' %%%
%%% in proc_, Normalized fft2 is applied in order to make A^{*} unitary. 
%%% SUBIM truncate the image at (1:SUBIM, 1:SUBIM); while proc_Nfft
%%% truncate the image by ratio. 1/SUMIM.
%%% 
%%% proc_TV has one more input 'str2', it combine two images together by 
%%% str1+i*str2.



%%% str1/ str2 :        location of the image.


function [IM,fft2_IM,fft2_IM_pois,SNR,mask]=proc_TV(str1,str2,os_rate,num_mask,add_phase,...
    sector_alpha,sector_beta,sector_mode,coer,SUBIM, rand_im)        

if nargin<11
    rand_im=0;
end

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

IM=IM1+1i*IM2 *0;   % construct str1+1i*str2.

IM=double(IM);

IM=IM/coer;        %%% coer to regulate poisson SNR. coer =400; SNR=3-4;
IM=min(IM,255);


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