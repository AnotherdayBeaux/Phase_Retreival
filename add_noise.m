
% IM is the image with random phase. 
% coer regulate the SNR, the noise added is poisson noise.

function [SNR,IM,fft2_IM,fft2_IM_pois]=add_noise(IM,coer,os_rate,num_mask,mask)

IM=IM/coer;

fft2_IM=Nos_fft(IM,os_rate,num_mask,mask);
  
fft2_IM_pois1=abs(fft2_IM).^2;
fft2_IM_pois=poissrnd(fft2_IM_pois1,size(fft2_IM));                     % poisson noise on sector
SNR=norm(abs(fft2_IM).^2,'fro')/norm(abs(fft2_IM).^2-fft2_IM_pois,'fro');


end