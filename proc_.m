%%% combination of 'add_noise'+'add_phase_mask' %%%
%%% in proc_, Normalized fft2 is applied in order to make A^{*} unitary. 
%%% SUBIM truncate the image at (1:SUBIM, 1:SUBIM); while proc_Nfft
%%% truncate the image by ratio. 1/SUMIM.



%%% str:        location of the image.


function [IM,fft2_IM,fft2_IM_pois,SNR,mask]=proc_(str,os_rate,num_mask,add_phase,...
    sector_alpha,sector_beta,sector_mode,coer,SUBIM, rand_im)        

if nargin<10
    rand_im=0;
end

IM=imread(str);
if length(size(IM))>2
    IM=rgb2gray(IM);   %%% convert the image to grayscale. 
end
IM=IM(1:SUBIM,1:(SUBIM));     %%% consider a sub image problem in order to save time.

ranS=rand(size(IM));ranS(find(ranS>0.2))=0;
IM=ranS*100;
%IM=double(IM)

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

if rand_im ~= 0
    IM=10*IM.*rand(size(IM)).*exp(1i*2*pi.*rand(size(IM)));
end



switch num_mask             % 0  0 mask / Id mask
                            % 1  1 mask case
                            % 2  2 mask case 
                            % 3  1 and 1/2 mask case
                            % 4  3 mask case, 1 id 2 random mask
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

fft2_IM=Nos_fft(IM,os_rate,num_mask,mask);
  
fft2_IM_pois1=abs(fft2_IM).^2;
fft2_IM_pois=poissrnd(fft2_IM_pois1,size(fft2_IM));                     % poisson noise on sector
SNR=norm(abs(fft2_IM).^2,'fro')/norm(abs(fft2_IM).^2-fft2_IM_pois,'fro');
end