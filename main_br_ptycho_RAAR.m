%%%%%% Blind random scan ptychography  no-noise case RAAR method





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
Tol=1e-7;

SUBIM=140;   % truncate image size cameraman default 256-by-256

bd_con=1;     % periodic boundary condition.
p_c=0;        % perturb center
type_im=1;
coer=1;


L=41; % measure the size of each small block (patch) L-by-L.



BETA=0.1:0.1:0.9;
PHASE_REG=0.1:0.1:0.7;
relative_DR_x=zeros(MaxIter,1);
resi_DR_y=zeros(MaxIter,1);
Success_Count= zeros(length(PHASE_REG),length(BETA));
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



Bigname='/Users/Beaux/Desktop/phase retrieval/code/16/br_ptycho_RAAR';
mkdir(Bigname);
x_line=1:20:SUBIM; y_line=1:20:SUBIM;
x_c_p=(1:20:SUBIM)'*ones(1,length(x_line));
y_c_p=ones(length(x_line),1)*(1:20:SUBIM);

x_c_p=x_c_p-unidrnd(7,size(x_c_p))+4;
y_c_p=y_c_p-unidrnd(7,size(y_c_p))+4;



    
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
              
count_phase_reg=1;
for phase_reg=PHASE_REG
    %%% lambda_0=Y+rand(size(Y))*4;
    lambda_0=10*(rand(size(Y))+1i*rand(size(Y)));
    
    est_phase=1/2*rand(l_patch,l_patch)-1/4;
    mask=exp(2i*pi*(phase_arg));
    mask_estimate= exp(2i*pi*(phase_arg+phase_reg*est_phase));
    Foldername=[Bigname,'/phase_arg=',num2str(phase_reg)];
    mkdir(Foldername)
    count_beta=1;

for beta = BETA
    %DR_Pois
    lambda_t=lambda_0;
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
    
    while  update_im < MaxIter_u_ip
    % phase retreival applied to image
            Plambda_t=B.*lambda_t./abs(lambda_t);
            Rlambda_t=2*Plambda_t-lambda_t;
    
            
            x_t=ptycho_rand_ifft(Rlambda_t,Na,Nb,os_rate,num_mask,mask_estimate,l_patch,x_c_p,y_c_p,p_c); 
            x__t=x_t./nor_ptycho;
            AARlambda=ptycho_rand_fft(x__t,os_rate,num_mask,mask_estimate,l_patch,x_c_p,y_c_p,bd_con);
            lambda_t=beta*(lambda_t+AARlambda-(2*beta-1)/beta * Plambda_t);
           
            %ee=norm(IM(:)'*x__t(:))/(IM(:)'*x__t(:));
            
            %rel_xx=norm(ee*x__t-IM,'fro')/norm(IM,'fro');
            
            res_y= norm(abs(Y)-abs(Rlambda_t),'fro')/norm_Y;
            fprintf('update_im=%d\n res_y=%s\n',update_im, res_y)
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
    
    %while  update_phase < MaxIter_u_ip

            

    %       Pfft_phase= B.*fft_phase./abs(fft_phase);
    %       Rfft_phase=2*Pfft_phase-fft_phase;
            
    %        nor_mask_tplus1=pr_phase_ifft(fft_phase,os_rate,num_mask,x__t,l_patch,x_c_p,y_c_p)./sqrt(nor_phase);
    %        mask_tplus1=nor_mask_tplus1./abs(nor_mask_tplus1);
    %        nor_mask_tplus1= sqrt(nor_phase).*mask_tplus1;
    %        fft_phase_1=pr_phase_fft(nor_mask_tplus1./sqrt(nor_phase),os_rate,num_mask,x__t,l_patch,x_c_p,y_c_p);
            
            
    %        fft_phase=beta*(fft_phase+fft_phase_1-(2*beta-1)/beta * Pfft_phase );
            
            
            
            
            %ee=norm(mask(:)'*mask_tplus1(:))/(mask(:)'*mask_tplus1(:));
            %rel_mask=norm(ee*mask_tplus1 - mask,'fro')/norm(mask,'fro');
     %       res_mask=norm(abs(fft_phase_1)-b,'fro')/norm_Y;
            %fprintf('u_p=%d\n rel_mask=%f\n res_mask=%f\n',update_phase,rel_mask,res_mask);
     %       fprintf('u_p=%d\n res_mask=%f\n',update_phase,res_mask);
     %       update_phase=update_phase+1;
            
    % end
    
    %%%%%%%
    %%%%%%% use pois likelihood ADMM 
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
    ratio=norm(mask_estimate,'fro')/norm(ones(size(mask_estimate)),'fro');
    
    
    
    ee=norm(IM(:)'*x__t(:))/(IM(:)'*x__t(:));
    rel_x=norm(ee*x__t-IM,'fro')/norm(IM,'fro');
    resi_y=norm(abs(AARlambda)-B,'fro')/norm(b,'fro');
    resi_DR_y(count_DR)= resi_y;
    relative_DR_x(count_DR)=rel_x;
    fprintf('count_DR=%d\n rel_x=%f\n', count_DR,rel_x);
    count_DR=count_DR+1;
    
    
    
    end
    
    
    foldername= [Foldername,'/beta=',num2str(beta)];
    mkdir(foldername)
    figure
    semilogy(1:count_DR-1,relative_DR_x(1:count_DR-1),'g')
    xlabel('iter')
    ylabel('Error')
    hold on
    semilogy(1:count_DR-1,resi_DR_y(1:count_DR-1),'r')
    legend('rel','resi')
    name=[foldername,'/error'];
    saveas(gcf,name,'png')
    close
 
    ee=norm(IM(:)'*x__t(:))/(IM(:)'*x__t(:));
    x_t_=ee*x__t;
    figure
    imshow(real(x_t_)/255)
    name=[foldername,'/rec_im'];
    saveas(gcf,name,'png')
    close
    if rel_x < 1e-4
        
        Success_Count(count_phase_reg,count_beta)=1;
    end
 count_beta=count_beta+1;   
    
end

count_phase_reg=count_phase_reg+1;
end






