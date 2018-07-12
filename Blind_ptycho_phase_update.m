%%%%% test the convergence of |A^x phase| = 0;

% i think ADMM will work globally, but raar will not




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

beta=0.2;
gamma=0.9;
relative_DR_x=zeros(MaxIter,1);
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



Bigname='/Users/Beaux/Desktop/phase retrieval/code/18/ptychography_2018_4_16';
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
              
 
% initialization of phase wrt phase constraint. 
est_phase=1/2*rand(l_patch,l_patch)-1/4;
mask=exp(2i*pi*(phase_arg));
mask_estimate= exp(2i*pi*(phase_arg+1*est_phase));
    
    
rel_x=1;
MaxIter_u_ip=100;   
rel_RAAR=zeros(MaxIter_u_ip,1);
rel_ADMM=zeros(MaxIter_u_ip,1);

    x__t= IM;
    update_phase_RAAR=1;
   % calculate nor_phase dependent on image x__t
    phase_=ones(l_patch,l_patch);
    fft_phase_=pr_phase_fft(phase_,os_rate,num_mask,x__t,l_patch,x_c_p,y_c_p); 
    nor_phase=pr_phase_ifft(fft_phase_,os_rate,num_mask,x__t,l_patch,x_c_p,y_c_p);


    new_mask = mask_estimate;
    nor_mask=new_mask.*sqrt(nor_phase);
    
    fft_phase_or=pr_phase_fft(nor_mask./sqrt(nor_phase),os_rate,num_mask,x__t,l_patch,x_c_p,y_c_p);
    fft_phase=fft_phase_or;
    
    while  update_phase_RAAR < MaxIter_u_ip

            

           Pfft_phase= B.*fft_phase./abs(fft_phase);
           Rfft_phase=2*Pfft_phase-fft_phase;
            
            nor_mask_tplus1=pr_phase_ifft(fft_phase,os_rate,num_mask,x__t,l_patch,x_c_p,y_c_p)./sqrt(nor_phase);
            mask_tplus1=nor_mask_tplus1./abs(nor_mask_tplus1);
            nor_mask_tplus1= sqrt(nor_phase).*mask_tplus1;
            fft_phase_1=pr_phase_fft(nor_mask_tplus1./sqrt(nor_phase),os_rate,num_mask,x__t,l_patch,x_c_p,y_c_p);
            
            
            fft_phase=beta*(fft_phase+fft_phase_1-(2*beta-1)/beta * Pfft_phase );
            
            
            
            
            ee=norm(mask(:)'*mask_tplus1(:))/(mask(:)'*mask_tplus1(:));
            rel_mask=norm(ee*mask_tplus1 - mask,'fro')/norm(mask,'fro');
            res_mask=norm(abs(fft_phase_1)-b,'fro')/norm_Y;
            fprintf('RAAR=%d\n rel_mask=%f\n res_mask=%f\n',update_phase_RAAR,rel_mask,res_mask);
            %fprintf('u_p=%d\n res_mask=%f\n',update_phase,res_mask);
            
            rel_RAAR(update_phase_RAAR)=rel_mask;
            update_phase_RAAR=update_phase_RAAR+1;
            
     end
    
    %%%%%%%
    update_phase_ADMM=1;
    %%%%%%% use pois likelihood ADMM 
    fft_phase= fft_phase_or;
    while  update_phase_ADMM < MaxIter_u_ip
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
            fprintf('ADMM=%d\n rel_mask=%f\n res_mask=%f\n',update_phase_ADMM,rel_mask,res_mask);
            rel_ADMM(update_phase_ADMM)=rel_mask;
            update_phase_ADMM=update_phase_ADMM+1;
            
    end
    
    figure
    semilogy(1:update_phase_RAAR-1, rel_RAAR(1:update_phase_RAAR-1),'r')
    hold on
    semilogy(1:update_phase_ADMM-1, rel_ADMM(1:update_phase_ADMM-1),'g')
    ylabel('RelEroor')
    ylabel('Iter')
    legend('RAAR','ADMM')
    title('phase update with precise image')
    
    
