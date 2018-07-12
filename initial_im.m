%%% combination of 'add_noise'+'add_phase_mask' %%%
%%% SUBIM truncate the image at (1:SUBIM, 1:SUBIM); while proc_Nfft
%%% truncate the image by ratio. 1/SUMIM.
%%% 
%%% proc_TV has one more input 'str2', it combine two images together by 
%%% str1+i*str2.

%%% type_im=0 real im
%%%        =1 str1+i*str2
%%%        =2 str1+i*str1



%%% str1  str2 :        location of the image.


function [IM,fft2_IM,fft2_IM_pois,SNR,mask]=initial_im(str1,str2,os_rate,num_mask,coer,SUBIM, type_im, RanS)        

if nargin==6
    type_im=2;
end

IM1=imread(str1);
if length(size(IM1))>2
    IM1=rgb2gray(IM1);   %%% convert the image to grayscale. 
end
IM1=IM1(1:SUBIM,1:(SUBIM));     %%% consider a sub image problem in order to save time.
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
end

IM=IM/coer;        %%% coer to regulate poisson SNR. coer =400; SNR=3-4;
%IM=ones(size(IM))*100/coer;

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
    case 4
        mask=ones([size(IM),3]);
        mask(:,:,1)=exp(2i*pi*rand(size(IM)));
        mask(:,:,2)=exp(2i*pi*rand(size(IM)));
end

fft2_IM=random_sample_fft(IM,os_rate,num_mask,mask,RanS);
  
fft2_IM_pois1=abs(fft2_IM).^2;
fft2_IM_pois=poissrnd(fft2_IM_pois1,size(fft2_IM));                     % poisson noise on sector
SNR=norm(abs(fft2_IM).^2,'fro')/norm(abs(fft2_IM).^2-fft2_IM_pois,'fro'); % snr of data.
end