%%%% ptychograph initialization. 

%%% subim = size of each small patch
%%% bd_con = 0 ; 0 boundary condition 
%%%        = 1 ; periodic boundary condtion
%%%        = 2 ; not consider boundary



function [enl_IM,fft2_IM,fft2_IM_pois,NSR,mask,nor_ptycho]=ini_ptycho(str1,str2,os_rate,...
    num_mask,coer,SUBIM, type_im,subim,overlap, bd_con)        

num_mask=1;



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
end

IM=IM/coer;        %%% coer to regulate poisson SNR. coer =400; SNR=3-4;
%IM=ones(size(IM))*100/coer;

switch num_mask             % 0  0 mask / Id mask
                            % 1  1 mask case
                            % 2  2 mask case 
                            % 3  1 and 1/2 mask case
    case 0
        mask = ones(subim,subim);
    case 1
        mask = exp(2i*pi*rand(subim,subim));
    case 2
        mask = exp(2i*pi*rand(subim,subim,2));
    case 3
        mask = ones(subim,subim,2);
        mask(:,:,1) = exp(2i*pi*rand(subim,subim));
    case 4
        mask=ones(subim,subim,3);
        mask(:,:,1)=exp(2i*pi*rand(subim,subim));
        mask(:,:,2)=exp(2i*pi*rand(subim,subim));
end




[Na,Nb]=size(IM);
patch=subim/2;
switch bd_con
    case 0    % 0 boundary
        enl_IM=zeros(subim+size(IM));
        enl_IM(patch+1:Na+patch, patch+1:Nb+patch)= IM;
    case 1    % periodic 
        Big_IM=kron(ones(3,3),IM);
        enl_IM=Big_IM(Na-patch+1:2*Na+patch,Nb-patch+1:2*Nb+patch);
        
    case 2
        enl_IM=IM;
end




IM_=ones(size(enl_IM));
p_fft_IM=ptycho_fft(IM_,os_rate,num_mask,mask,subim,overlap);
IM__=ptycho_ifft(p_fft_IM,os_rate,num_mask,mask,subim,overlap);
nor_ptycho=sqrt(IM__);


fft2_IM=ptycho_fft(enl_IM,os_rate,num_mask,mask,subim,overlap);
  
fft2_IM_pois1=abs(fft2_IM).^2;
fft2_IM_pois=poissrnd(fft2_IM_pois1,size(fft2_IM));                     % poisson noise on sector
NSR=norm(abs(fft2_IM)-sqrt(fft2_IM_pois),'fro')/norm(abs(fft2_IM),'fro'); 
% NSR of data. agree with Chen. 
end


