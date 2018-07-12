%%% add phase to original image (trancated, 1/SUBIM of IM), and generate
%%% masks. 

%%% str = location of image.
%%% num_mask = 0  0 mask / Id mask
             % 1  1 mask case
             % 2  2 mask case 
             % 3  1 and 1/2 mask case
%%% add_phase = 1 % add phase
%%%           = 0 % real image
%%% sector_mode = look below.
%%% SUBIM means pick 1/SUBIM of IM in both x-axis and y-axis. 

%%% size(mask)= [Na,Nb,1] when num_mask = 0,1
%%%           = [Na,Nb,2] when num_mask = 2,3









function [IM,mask]=add_phase_mask(str,num_mask,add_phase,sector_alpha,sector_beta,sector_mode,SUBIM)        


IM=imread(str);
if length(size(IM))>2
    IM=rgb2gray(IM);   %%% convert the image to grayscale. 
end
[sub_x,sub_y]=size(IM);

IM=IM(1:floor(sub_x/SUBIM),1:floor(sub_y/SUBIM));     %%% consider a sub image problem in order to save time.
IM=double(IM);        


if add_phase == 1
    switch sector_mode  %%%  sector_mode 0 general sector constrains; [-alpha pi, beta pi]%
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

end