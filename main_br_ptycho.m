%%%%%% Blind  ptychography  no-noise case





Bigname='/Users/Beaux/Desktop/phase retrieval'   ; % location of image

str1=[Bigname,'/image_lib/Cameraman.bmp'];   %% consider Im= Cameraman+i* barbera
str2=[Bigname,'/image_lib/Cameraman.bmp'];

os_rate=2;  %% oversampling ratio

            
 
      
num_mask=1; % 0  0 mask / Id mask
            % 1  1 mask case
            % 2  2 mask case 
            % 3  1 and 1/2 mask case
            
                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                 %%%
%%%  sector constraint parameter    %%%
%%%  x_k \in [-alpha \pi, beta \pi] %%%
%%%                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  sector_mode 0 general sector constrains; %
%%%              1 real  positivity;          %
%%%              2 complex positivity;        % 
%%%              3 non constrains;  

sector_alpha= 0.1;
sector_beta=0.1;
sector_mode=3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% add phase to the original image %%%
%%%                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The phase added should be consistent
%%% with the sector constrain.     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

add_phase=1;% 0 real image
            % 1 add phase consistent with sector constrain
            
MaxIter=400;
Tol=1e-14;

SUBIM=256;   % truncate image size cameraman default 256-by-256

bd_con=1;     % periodic boundary condition.
p_c=0;        % perturb center
type_im=1;
coer=1;


L=65; % measure the size of each small block (patch) L-by-L.

r=0.3; % phase mask constraint r< 1/2 always.
gamma=1;
alpha=0.0001;
relative_DR_y=zeros(MaxIter,1);
resi_DR_y=zeros(MaxIter,1);
len_l=length(L);


%%% AMDD method 
%%% L(x,y,lambda) = sum{loglikelyhood + lagrange multiplier + p/2 *augmented LM} 
%%% sub to y=Ax
%%% 1st  y_k+1 = argmin L(x_k,y,lambda_k)
%%% 2nd  x_k+1 = argmin L(x,y_k+1,lambda_k)
%%% 3rd  lambda_k+1=lambda_k+p(y_k+1-Ax_k+1)


%%%%%%%%%%%%%%%%%%%%%%
%%%                %%%
%%% Initialization %%%
%%%                %%%
%%%%%%%%%%%%%%%%%%%%%%



Bigname='/Users/Beaux/Desktop/phase retrieval/code/16/ptychography_2018_4_16';
mkdir(Bigname);
x_line=1:30:256; y_line=1:30:256;
x_c_p=(1:30:256)'*ones(1,length(x_line));
y_c_p=ones(length(x_line),1)*(1:30:256);

%x_c_p=x_c_p-unidrnd(7,size(x_c_p))+4;
%y_c_p=y_c_p-unidrnd(7,size(y_c_p))+4;



    
l_patch=L;
%ini_rand_ptycho(str1,str2,os_rate,...
%   num_mask,coer,SUBIM, type_im,l_patch,x_c_p, y_c_p,bd_con) 
[IM,Y,Z,NSR,phase_arg,nor_ptycho,nor_phase]=ini_rand_ptycho(str1,str2,os_rate,num_mask,coer,SUBIM,...
    type_im,l_patch,x_c_p,y_c_p,bd_con);
 
 
[Na,Nb]=size(IM);
norm_IM=norm(IM,'fro');

[Y_Na,Y_Nb]=size(Y);
    
Z=abs(Y).^2;
b=abs(Y);
norm_Y=norm(b,'fro');
B=sqrt(abs(Z));  %%% noise case
              
 

    %%% lambda_0=Y+rand(size(Y))*4;
    lambda_0=10*(rand(size(Y))+1i*rand(size(Y)));

    %DR_Pois
    lambda_t=lambda_0;
    est_phase=1/2*rand(l_patch,l_patch)-1/4;
    mask=exp(2i*pi*(phase_arg));
    mask_estimate= exp(2i*pi*(phase_arg+1*est_phase));
    count_DR=1;   
    rel_x=1;
    
    
    while rel_x > Tol && count_DR< MaxIter
    
    MaxIter_u_ip=(count_DR < 15)*40+(count_DR >= 15 && count_DR < 30 )* 10 + (count_DR>= 45) *5;
    
    update_im=1;
    %%% calculate nor_factors depenent on image
    IM_=ones(size(IM));
    [Na,Nb]=size(IM_);
    p_fft_IM=ptycho_rand_fft(IM_,os_rate,num_mask,mask_estimate,l_patch,x_c_p,y_c_p,bd_con);
    nor_ptycho=ptycho_rand_ifft(p_fft_IM,Na,Nb,os_rate,num_mask,mask_estimate,l_patch,x_c_p,y_c_p,p_c);
    %%% end calculate nor_factors
    
    residual_x=zeros(1,MaxIter_u_ip);
    while  update_im < MaxIter_u_ip
    % phase retreival applied to image
            x_t=ptycho_rand_ifft(lambda_t,Na,Nb,os_rate,num_mask,mask_estimate,l_patch,x_c_p,y_c_p,p_c); 
            x__t=x_t./nor_ptycho;
            AAlambda=ptycho_rand_fft(x__t,os_rate,num_mask,mask_estimate,l_patch,x_c_p,y_c_p,bd_con);
            y_tplus1=AAlambda;

            ylambda_k=2*y_tplus1-lambda_t;
            Q=abs(ylambda_k);
            % POISSON
            rot_z=(Q/gamma+sqrt(Q.^2/gamma^2+8*(2+1/gamma)*(Z)))/(4+2/gamma);
            
            % GAUSSIAN
            %rot_z=(sqrt(Z)+1/gamma*Q)/(1+1/gamma);
            
            z_tplus1=rot_z.ylambda_k./Q;
         
            lambda_t=lambda_t+z_tplus1-y_tplus1;
            ee=norm(IM(:)'*x__t(:))/(IM(:)'*x__t(:));

            rel_xx=norm(ee*x__t-IM,'fro')/norm(IM,'fro');
            res_x= norm(abs(Y)-abs(y_tplus1),'fro')/norm_Y;
            fprintf('u_i=%d\n rel_x=%f\n res_x=%f\n', update_im,rel_xx,res_x);
            residual_x(update_im)=res_x;
            update_im=update_im+1;
     
    end
    %%%%%%%%%%
    %%%%%%%%%%
    %figure
    %semilogy(1:update_im-1,residual_x(1,1:update_im-1),'g')
    %xlabel('iter')
    %ylabel('resError')
    %%%%%%%%%%%%%%%%%%%
    
   % phase retreival applied to mask
    update_phase=1;
   % calculate nor_phase dependent on image
    phase_=ones(l_patch,l_patch);
    fft_phase_=pr_phase_fft(phase_,os_rate,num_mask,x__t,l_patch,x_c_p,y_c_p); 
    nor_phase=pr_phase_ifft(fft_phase_,os_rate,num_mask,x__t,l_patch,x_c_p,y_c_p);


    new_mask = mask_estimate;
    nor_mask=new_mask.*sqrt(nor_phase);
    
    fft_phase=pr_phase_fft(nor_mask./sqrt(nor_phase),os_rate,num_mask,x__t,l_patch,x_c_p,y_c_p);
    
    while  update_phase < MaxIter_u_ip
            nor_mask_tplus1=pr_phase_ifft(fft_phase,os_rate,num_mask,x__t,l_patch,x_c_p,y_c_p)./sqrt(nor_phase);
            mask_tplus1=nor_mask_tplus1./abs(nor_mask_tplus1);
            nor_mask_tplus1= sqrt(nor_phase).*mask_tplus1;
            fft_phase_1=pr_phase_fft(nor_mask_tplus1./sqrt(nor_phase),os_rate,num_mask,x__t,l_patch,x_c_p,y_c_p);
            
            Q_phase=2*fft_phase_1-fft_phase;
            Q=abs(Q_phase);
            rot_z=(Q/gamma+sqrt(Q.^2/gamma^2+8*(2+1/gamma)*(Z)))/(4+2/gamma);
            z_mask_tplus1 = rot_z.* Q_phase./Q;
            
            
            fft_phase=fft_phase+ z_mask_tplus1   -fft_phase_1;
            
            
            ee=norm(mask(:)'*mask_tplus1(:))/(mask(:)'*mask_tplus1(:));
            rel_mask=norm(ee*mask_tplus1 - mask,'fro')/norm(mask,'fro');
            res_mask=norm(abs(fft_phase_1)-b,'fro')/norm_Y;
            fprintf('u_p=%d\n rel_mask=%f\n res_mask=%f\n',update_phase,rel_mask,res_mask);
            update_phase=update_phase+1;
            
    end
    mask_estimate=nor_mask_tplus1./sqrt(nor_phase);
    
    
    
    ee=norm(IM(:)'*x__t(:))/(IM(:)'*x__t(:));
    eee=norm(Y(:)'*lambda_t(:))/(Y(:)'*lambda_t(:));
    rel_y=norm(eee*lambda_t-Y,'fro')/norm(Y,'fro');
    rel_x=norm(ee*x__t-IM,'fro')/norm(IM,'fro');
    
    resi_DR_y(1,count_DR)= res_x;
    relative_DR_y(1,count_DR)=rel_x;
    fprintf('count_DR=%d\n rel_x=%f\n', count_DR,rel_x);
    count_DR=count_DR+1;
    
    
    
    end
    
    
    
    
    figure
    semilogy(1:count_DR-1,relative_DR_y(1,1:count_DR-1),'g')
    xlabel('iter')
    ylabel('Error')
    hold on
    semilogy(1:count_DR-1,resi_DR_y(1,1:count_DR-1),'r')
    legend('rel','resi')
    
 
    ee=norm(IM(:)'*x__t(:))/(IM(:)'*x__t(:));
    x_t_=ee*x__t;
    figure
    imshow(real(x_t_)/255)
