%%%% Blind  ptychography  no-noise case

%%%% |phi(n)-phi_0| < \delta * \pi; phi_0 is uniform for all element in the
%%%% mask

%%%% The main task of this program is to
%%%%     1. explore the effect of different mask a) iid mask 
%                                                b) correlated mask
%                                                c) frensel mask




Bigname='/Users/Beaux/Desktop/phase retrieval'   ; % location of image

str1=[Bigname,'/image_lib/Cameraman.bmp'];   %% consider Im= Cameraman+i* barbera
str2=[Bigname,'/image_lib/Barbara256.png'];

os_rate=2;  %% oversampling ratio

            
      
num_mask=1; % 0  0 mask / Id mask 
            % 1  1 mask case
            % 2  2 mask case 
            % 3  1 and 1/2 mask case
            
                         

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% add phase to the original image %%%
%%%                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The phase added should be consistent
%%% with the sector constrain.     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         
MaxIter=200;
Tol=1e-13;

SUBIM=256;   % truncate image size cameraman default 256-by-256

bd_con=1;     % periodic boundary condition.
p_c=0;        % perturb center
type_im=2;    % str1+i*str2 
mask_type=1;  % 2 fresnel mask
              % 1 correlated mask
              % 0 iid mask
              
Likelihood=0;  % poisson ==0
               % gaussian ==1


              
coer=1;    % regulate noise poisson type




gamma=1;
relative_DR_y=zeros(MaxIter,5);
resi_DR_y=zeros(MaxIter,5);
o_l = 0.5;%0.2:0.1:0.8 ; % overlaping
cim_diff=30; %20:15:120;
perturb = 3;    % set the max magnitude of either horizontal or vertical perturb
rank_pert = 100;% =1 rank_one perturb
  % delta < 1/2 regulates the MPC.
delta=1/2;
delta_cor=1;
cor_m=1;
Count_DR=zeros(5,1);
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
% cim_diff = center_image_difference  regulates the number of small patch.
% l_patch  =side length of patch      regulates the size of each small
% patch.
    
    
     l_patch_pre= (cim_diff+1)/(1-o_l);
     l_patch = floor(l_patch_pre/2) *2+1;
     q=floor(SUBIM/cim_diff+1);
     num_subim=q^2;  
 
     beta_1=l_patch/2; beta_2 = l_patch/2; rho =0.5;
     
    %%%%%% change the rank condition
     
     
     

    x_line=1:cim_diff:SUBIM; y_line=1:cim_diff:SUBIM;
    x_c_p=(1:cim_diff:SUBIM)'*ones(1,length(x_line));
    y_c_p=ones(length(x_line),1)*(1:cim_diff:SUBIM);
    
    

    if rank_pert==1
    % rank one; set specific perturb pattern
    x_c_p=x_c_p-(unidrnd(perturb*2+1,size(x_c_p,1),1)-perturb-1)*(ones(1,size(x_c_p,2)));
    y_c_p=y_c_p-(ones(size(y_c_p,2),1))*(unidrnd(perturb*2+1,1,size(y_c_p,2))-perturb-1);
    % full rank perturbation
    else
    x_c_p=x_c_p-unidrnd(perturb*2+1,size(x_c_p))+perturb+1;
    y_c_p=y_c_p-unidrnd(perturb*2+1,size(y_c_p))+perturb+1;

    end
    Bigname='/Users/Beaux/Desktop/phase retrieval/code/16/diff_cor';
    mkdir(Bigname)
    
    
    cl_count=1;
    
    
    %%%% phase_arg fixed for all correlated and iid mask case %%%%%
    phase_arg_fix=rand(l_patch,l_patch);
    %%%%
    
    
    phase_arg_input=phase_arg_fix;   %%% make sure phase_arg is an iid mask
    
    
    
    lambda_0=10*(rand(size(Y))+1i*rand(size(Y))); % random initialize  image
    est_phase=delta*(rand(l_patch,l_patch)-1/2);  % perturb mask initialization
    
    
    
    
    
    for c_l_ratio = 0.4:0.3:0.7
        
        c_l=floor(l_patch*c_l_ratio);
    
    
    
    [IM,Y,Z,NSR,phase_arg,nor_ptycho,nor_phase]=ini_rand_ptycho_argfixed(str1,str2,os_rate,num_mask,coer,SUBIM,...
    type_im,l_patch,x_c_p,y_c_p,bd_con,mask_type, phase_arg_input,c_l, beta_1,beta_2, rho);
 
 
    [Na,Nb]=size(IM);
    norm_IM=norm(IM,'fro');

    [Y_Na,Y_Nb]=size(Y);
    
    b=abs(Y);
    norm_Y=norm(b,'fro');
    Z=abs(Y).^2;% no-noise 
    

    
    
   
    mask=exp(2i*pi*(phase_arg));
    %DR_Pois
    
    
    
    %%%%%%% pick the correct likelihood function
    if Likelihood == 0
        Rot_z = @(Q) (Q/gamma+sqrt(Q.^2/gamma^2+8*(2+1/gamma)*(Z)))/(4+2/gamma); 
    elseif Likelihood ==1
        Rot_z = @(Q) (sqrt(Z)+1/gamma*Q)/(1+1/gamma);
    end
    
    lambda_t=lambda_0;
    
    %%%%% correlated perturbation.
    
    
    %mask_estimate= exp(2i*pi*phase_arg).*l_phase_drift_mask.*  exp(1/2*2i*pi*(rand(l_patch,l_patch)-1/2));
    %mask_estimate= exp(2i*pi*phase_arg).*l_phase_drift_mask.*  Cor_Mask(m,c_l,1/2);
    
    mask_estimate= exp(2i*pi*phase_arg).* exp(2i*pi*est_phase) ;
    %%%%%
   
    
    mask_estimate_ini=mask_estimate;
    rel_x=1;
    count_DR=1;
while rel_x > Tol && count_DR< MaxIter
     MaxIter_u_ip = 60;
    %MaxIter_u_ip=(count_DR==1)*60+(count_DR < 20 && count_DR >=2)*40+...
     %   (count_DR >= 20 && count_DR < 40 )* 20 +...
     %   (count_DR>= 40 && count_DR<60) *10+...
     %   (count_DR >= 60 ) *5;
    
    update_im=1;
    %%% calculate nor_factors depenent on image
    IM_=ones(size(IM));
    [Na,Nb]=size(IM_);
    p_fft_IM=ptycho_rand_fft(IM_,os_rate,num_mask,mask_estimate,l_patch,x_c_p,y_c_p,bd_con);
    nor_ptycho=ptycho_rand_ifft(p_fft_IM,Na,Nb,os_rate,num_mask,mask_estimate,l_patch,x_c_p,y_c_p,p_c);
    %%% end calculate nor_factors
    
    residual_x=zeros(1,MaxIter_u_ip+3);
    resi_diff_im=1;
    x__tpre=0;
    while  update_im < 5 || ( update_im < MaxIter_u_ip && resi_diff_im >1e-4)
    % phase retreival applied to image
            x_t=ptycho_rand_ifft(lambda_t,Na,Nb,os_rate,num_mask,mask_estimate,l_patch,x_c_p,y_c_p,p_c); 
            x__t=x_t./nor_ptycho;
            AAlambda=ptycho_rand_fft(x__t,os_rate,num_mask,mask_estimate,l_patch,x_c_p,y_c_p,bd_con);
            y_tplus1=AAlambda;

            ylambda_k=2*y_tplus1-lambda_t;
            Q=abs(ylambda_k);
            rot_z=Rot_z(Q);
            
            z_tplus1=rot_z.*ylambda_k./Q;
           
            lambda_t=lambda_t+z_tplus1-y_tplus1;
  
            % end test
            %ee=norm(IM(:)'*x__t(:))/(IM(:)'*x__t(:));
            
            %rel_xx=norm(ee*x__t-IM,'fro')/norm(IM,'fro');
            res_x= norm(abs(Y)-abs(y_tplus1),'fro')/norm_Y;
            %fprintf('u_i=%d\n rel_x=%f\n res_x=%f\n', update_im,rel_xx,res_x);
            %fprintf('u_i=%d\n rel_q=%f rel_x=%f res_x=%f\n', update_im, rel_q,rel_x,res_x);
            residual_x(update_im+3)=res_x;
            resi_diff_im = 1/3*norm(residual_x(update_im+1:update_im+3)-residual_x(update_im:update_im+2),1)/res_x;
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
    
    residual_xx=zeros(1,MaxIter_u_ip+2);
    resi_diff_phase=1;
    while  update_phase < 5 || ( update_phase < MaxIter_u_ip && resi_diff_phase > 1e-4)
            nor_mask_tplus1=pr_phase_ifft(fft_phase,os_rate,num_mask,x__t,l_patch,x_c_p,y_c_p)./sqrt(nor_phase);
            mask_tplus1=nor_mask_tplus1./abs(nor_mask_tplus1);
            nor_mask_tplus1= sqrt(nor_phase).*mask_tplus1;
            fft_phase_1=pr_phase_fft(nor_mask_tplus1./sqrt(nor_phase),os_rate,num_mask,x__t,l_patch,x_c_p,y_c_p);
            
            Q_phase=2*fft_phase_1-fft_phase;
            Q=abs(Q_phase);
            rot_z= Rot_z(Q);
            
            z_mask_tplus1 = rot_z.* Q_phase./Q;
            
            
            fft_phase=fft_phase+ z_mask_tplus1   -fft_phase_1;
            
            
            %ee=norm(mask(:)'*mask_tplus1(:))/(mask(:)'*mask_tplus1(:));
            %rel_mask=norm(ee*mask_tplus1 - mask,'fro')/norm(mask,'fro');
            res_mask=norm(abs(fft_phase_1)-b,'fro')/norm_Y;
            residual_xx(update_phase+3) = res_mask;
            resi_diff_phase = 1/3*norm(residual_xx(update_phase+1:update_phase+3)-residual_xx(update_phase:update_phase+2),1)/res_mask;
            
            %fprintf('u_p=%d\n rel_mask=%f\n res_mask=%f\n',update_phase,rel_mask,res_mask);
            update_phase=update_phase+1;
            
    end
    mask_estimate=nor_mask_tplus1./sqrt(nor_phase);
    
    
    %%%% MPC enforcement.
    %mask_arg=angle(mask_estimate./mask);       % center
    %center_arg=sum(sum(mask_arg))/m/m;
    %proj_mask_arg=max(mask_arg-center_arg, -delta*pi);
    %proj_mask_arg=min(proj_mask_arg, delta*pi);
    %mask_estimate = exp(1i*proj_mask_arg);
    %%%% MPC enforcement.
    
    
    ee=norm(IM(:)'*x__t(:))/(IM(:)'*x__t(:));
    rel_x=norm(ee*x__t-IM,'fro')/norm(IM,'fro');
    
    resi_DR_y(count_DR,cl_count)= res_x;
    relative_DR_y(count_DR,cl_count)=rel_x;
    
    fprintf('count_DR=%d\n rel_x=%4.3s res_x=%s\n',...
                count_DR,rel_x, res_x);
    count_DR=count_DR+1;
    
    
   
end 

Count_DR(cl_count)=count_DR;



foldername=[Bigname,'/masktype=',num2str(mask_type),',c_l=',num2str(c_l_ratio)];
mkdir(foldername)

figure
semilogy(1:count_DR-1,relative_DR_y(1:count_DR-1,cl_count),'g')
hold on 
semilogy(1:count_DR-1,resi_DR_y(1:count_DR-1,cl_count),'g:')
hold on 
xlabel('Epoch')
ylabel(' ')
name=[foldername,'/convergence'];
legend('RE','Residual')
saveas(gcf,name,'png')
close
    
ee=norm(IM(:)'*x__t(:))/(IM(:)'*x__t(:));
x_t_=ee*x__t;
phase_diff=x_t_./IM;
ang_diff=angle(phase_diff);
figure
imshow(real(x_t_)/255)
title('real part of x_t')
name=[foldername,'/recon_RE'];
saveas(gcf,name,'png')
close
    
    
    figure
    imshow(imag(x_t_)/255)
    title('imag part of x_t')
    name=[foldername,'/recon_IM'];
    saveas(gcf,name,'png')
    close
    
    figure
    colormap('hot')
    imagesc(ang_diff)
    colorbar
    title('colormapIM')
    name=[foldername,'/ColormapIMDifference'];
    saveas(gcf,name,'png')
    close
    
    figure
    ee_mask=norm(mask_estimate_ini(:)'*mask(:))/(mask_estimate_ini(:)'*mask(:));
    mask_diff=ee_mask'*mask_estimate_ini./mask;
    mask_arg_diff=angle(mask_diff);
    colormap('hot')
    imagesc(mask_arg_diff)
    colorbar
    title('')
    name=[foldername,'/ColormapIni-MaskDifference'];
    h=gca;
    set(h,'FontSize',14)
    saveas(gcf,name,'png')
    close
    
    figure
    ee_mask=norm(mask_estimate(:)'*mask(:))/(mask_estimate(:)'*mask(:));
    mask_diff=ee_mask'*mask_estimate./mask;
    mask_arg_diff=angle(mask_diff);
    colormap('hot')
    imagesc(mask_arg_diff)
    %set(gcf,'Position',[500,300,256,256])
    colorbar
    title('colormapMASK')
    name=[foldername,'/ColormapMaskDifference'];
    saveas(gcf,name,'png')
    close
    
    %%%%% record mask angel distribution 0) iid 1)correlated 2) fresnel
    figure
    colormap('hot')
    imagesc(angle(mask))
    %set(gcf,'Position',[500,300,472,400])
    colorbar
    title('mask angle')
    name=[foldername,'/mask_angle'];
    h=gca;
    set(h,'FontSize',14);
    saveas(gcf,name,'png')
    close
    
    cl_count=cl_count+1;
    end
    
 c_l=l_patch;   
 for mask_type = 0:2
        

    
    [IM,Y,Z,NSR,phase_arg,nor_ptycho,nor_phase]=ini_rand_ptycho_argfixed(str1,str2,os_rate,num_mask,coer,SUBIM,...
    type_im,l_patch,x_c_p,y_c_p,bd_con,mask_type,phase_arg_input,c_l, beta_1,beta_2, rho);
 
 
    [Na,Nb]=size(IM);
    norm_IM=norm(IM,'fro');

    [Y_Na,Y_Nb]=size(Y);
    
    b=abs(Y);
    norm_Y=norm(b,'fro');
    Z=abs(Y).^2;% no-noise 
    
    
    
    
    mask=exp(2i*pi*(phase_arg));
    %DR_Pois
    
    
    
    %%%%%%% pick the correct likelihood function
    if Likelihood == 0
        Rot_z = @(Q) (Q/gamma+sqrt(Q.^2/gamma^2+8*(2+1/gamma)*(Z)))/(4+2/gamma); 
    elseif Likelihood ==1
        Rot_z = @(Q) (sqrt(Z)+1/gamma*Q)/(1+1/gamma);
    end
    
    lambda_t=lambda_0;
    
    %%%%% correlated perturbation.
    
    
    %mask_estimate= exp(2i*pi*phase_arg).*l_phase_drift_mask.*  exp(1/2*2i*pi*(rand(l_patch,l_patch)-1/2));
    %mask_estimate= exp(2i*pi*phase_arg).*l_phase_drift_mask.*  Cor_Mask(m,c_l,1/2);
    
    mask_estimate= exp(2i*pi*phase_arg).* exp(2i*pi*est_phase) ;
    %%%%%
   
    
    
    %%%% MPC enforcement.
    %mask_arg=angle(mask_estimate./mask);       % center
    %center_arg=sum(sum(mask_arg))/m/m;
    %proj_mask_arg=max(mask_arg-center_arg, -delta*pi);
    %proj_mask_arg=min(proj_mask_arg, delta*pi);
    %mask_estimate = exp(1i*proj_mask_arg);
    %%%% MPC enforcement.
    
   
    count_DR=1;   
    %rel_x=1;
    
    
    
    
    
    mask_estimate_ini=mask_estimate;
    rel_x=1;
while rel_x > Tol && count_DR< MaxIter
     MaxIter_u_ip = 60;
    %MaxIter_u_ip=(count_DR==1)*60+(count_DR < 20 && count_DR >=2)*40+...
     %   (count_DR >= 20 && count_DR < 40 )* 20 +...
     %   (count_DR>= 40 && count_DR<60) *10+...
     %   (count_DR >= 60 ) *5;
    
    update_im=1;
    %%% calculate nor_factors depenent on image
    IM_=ones(size(IM));
    [Na,Nb]=size(IM_);
    p_fft_IM=ptycho_rand_fft(IM_,os_rate,num_mask,mask_estimate,l_patch,x_c_p,y_c_p,bd_con);
    nor_ptycho=ptycho_rand_ifft(p_fft_IM,Na,Nb,os_rate,num_mask,mask_estimate,l_patch,x_c_p,y_c_p,p_c);
    %%% end calculate nor_factors
    
    residual_x=zeros(1,MaxIter_u_ip+3);
    resi_diff_im=1;
    x__tpre=0;
    while  update_im < 5 || ( update_im < MaxIter_u_ip && resi_diff_im >1e-4)
    % phase retreival applied to image
            x_t=ptycho_rand_ifft(lambda_t,Na,Nb,os_rate,num_mask,mask_estimate,l_patch,x_c_p,y_c_p,p_c); 
            x__t=x_t./nor_ptycho;
            AAlambda=ptycho_rand_fft(x__t,os_rate,num_mask,mask_estimate,l_patch,x_c_p,y_c_p,bd_con);
            y_tplus1=AAlambda;

            ylambda_k=2*y_tplus1-lambda_t;
            Q=abs(ylambda_k);
            rot_z=Rot_z(Q);
            
            z_tplus1=rot_z.*ylambda_k./Q;
           
            lambda_t=lambda_t+z_tplus1-y_tplus1;
  
            % end test
            %ee=norm(IM(:)'*x__t(:))/(IM(:)'*x__t(:));
            
            %rel_xx=norm(ee*x__t-IM,'fro')/norm(IM,'fro');
            res_x= norm(abs(Y)-abs(y_tplus1),'fro')/norm_Y;
            %fprintf('u_i=%d\n rel_x=%f\n res_x=%f\n', update_im,rel_xx,res_x);
            %fprintf('u_i=%d\n rel_q=%f rel_x=%f res_x=%f\n', update_im, rel_q,rel_x,res_x);
            residual_x(update_im+3)=res_x;
            resi_diff_im = 1/3*norm(residual_x(update_im+1:update_im+3)-residual_x(update_im:update_im+2),1)/res_x;
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
    
    residual_xx=zeros(1,MaxIter_u_ip+2);
    resi_diff_phase=1;
    while  update_phase < 5 || ( update_phase < MaxIter_u_ip && resi_diff_phase > 1e-4)
            nor_mask_tplus1=pr_phase_ifft(fft_phase,os_rate,num_mask,x__t,l_patch,x_c_p,y_c_p)./sqrt(nor_phase);
            mask_tplus1=nor_mask_tplus1./abs(nor_mask_tplus1);
            nor_mask_tplus1= sqrt(nor_phase).*mask_tplus1;
            fft_phase_1=pr_phase_fft(nor_mask_tplus1./sqrt(nor_phase),os_rate,num_mask,x__t,l_patch,x_c_p,y_c_p);
            
            Q_phase=2*fft_phase_1-fft_phase;
            Q=abs(Q_phase);
            rot_z= Rot_z(Q);
            
            z_mask_tplus1 = rot_z.* Q_phase./Q;
            
            
            fft_phase=fft_phase+ z_mask_tplus1   -fft_phase_1;
            
            
            %ee=norm(mask(:)'*mask_tplus1(:))/(mask(:)'*mask_tplus1(:));
            %rel_mask=norm(ee*mask_tplus1 - mask,'fro')/norm(mask,'fro');
            res_mask=norm(abs(fft_phase_1)-b,'fro')/norm_Y;
            residual_xx(update_phase+3) = res_mask;
            resi_diff_phase = 1/3*norm(residual_xx(update_phase+1:update_phase+3)-residual_xx(update_phase:update_phase+2),1)/res_mask;
            
            %fprintf('u_p=%d\n rel_mask=%f\n res_mask=%f\n',update_phase,rel_mask,res_mask);
            update_phase=update_phase+1;
            
    end
    mask_estimate=nor_mask_tplus1./sqrt(nor_phase);
    
    
    %%%% MPC enforcement.
    %mask_arg=angle(mask_estimate./mask);       % center
    %center_arg=sum(sum(mask_arg))/m/m;
    %proj_mask_arg=max(mask_arg-center_arg, -delta*pi);
    %proj_mask_arg=min(proj_mask_arg, delta*pi);
    %mask_estimate = exp(1i*proj_mask_arg);
    %%%% MPC enforcement.
    
    
    ee=norm(IM(:)'*x__t(:))/(IM(:)'*x__t(:));
    rel_x=norm(ee*x__t-IM,'fro')/norm(IM,'fro');
    
    resi_DR_y(count_DR,mask_type+3)= res_x;
    relative_DR_y(count_DR,mask_type+3)=rel_x;
    
    fprintf('count_DR=%d\n rel_x=%4.3s res_x=%s\n',...
                count_DR,rel_x, res_x);
    count_DR=count_DR+1;
    
    
   
end 

Count_DR(mask_type+3)=count_DR;



foldername=[Bigname,'/masktype=',num2str(mask_type)];
mkdir(foldername)

figure
semilogy(1:count_DR-1,relative_DR_y(1:count_DR-1,mask_type+3),'g')
hold on 
semilogy(1:count_DR-1,resi_DR_y(1:count_DR-1,mask_type+3),'g:')
hold on 
xlabel('Epoch')
ylabel(' ')
name=[foldername,'/convergence'];
legend('RE','Residual')
saveas(gcf,name,'png')
close
    
ee=norm(IM(:)'*x__t(:))/(IM(:)'*x__t(:));
x_t_=ee*x__t;
phase_diff=x_t_./IM;
ang_diff=angle(phase_diff);
figure
imshow(real(x_t_)/255)
title('real part of x_t')
name=[foldername,'/recon_RE'];
saveas(gcf,name,'png')
close
    
    
    figure
    imshow(imag(x_t_)/255)
    title('imag part of x_t')
    name=[foldername,'/recon_IM'];
    saveas(gcf,name,'png')
    close
    
    figure
    colormap('hot')
    imagesc(ang_diff)
    colorbar
    title('colormapIM')
    name=[foldername,'/ColormapIMDifference'];
    saveas(gcf,name,'png')
    close
    
    figure
    ee_mask=norm(mask_estimate_ini(:)'*mask(:))/(mask_estimate_ini(:)'*mask(:));
    mask_diff=ee_mask'*mask_estimate_ini./mask;
    mask_arg_diff=angle(mask_diff);
    colormap('hot')
    imagesc(mask_arg_diff)
    colorbar
    title('colormapIni-MASK')
    name=[foldername,'/ColormapIni-MaskDifference'];
    saveas(gcf,name,'png')
    close
    
    figure
    ee_mask=norm(mask_estimate(:)'*mask(:))/(mask_estimate(:)'*mask(:));
    mask_diff=ee_mask'*mask_estimate./mask;
    mask_arg_diff=angle(mask_diff);
    colormap('hot')
    imagesc(mask_arg_diff)
    colorbar
    title('colormapMASK')
    name=[foldername,'/ColormapMaskDifference'];
    saveas(gcf,name,'png')
    close
    
    %%%%% record mask angel distribution 0) iid 1)correlated 2) fresnel
    figure
    colormap('hot')
    imagesc(angle(mask))
    %set(gcf,'Position',[500,300,472,400])
    colorbar
    title('mask angle')
    name=[foldername,'/mask_angle'];
    h=gca;
    set(h,'FontSize',14)
    saveas(gcf,name,'png')
    close
    
    
end   
    

figure


semilogy(1,relative_DR_y(1,1),'r^:','LineWidth',2,...
                       'MarkerEdgeColor','r',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',10)
hold on
semilogy(2,relative_DR_y(2,2),'ro:','LineWidth',2,...
                       'MarkerEdgeColor','r',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',10)
hold on
semilogy(3,relative_DR_y(3,4),'r*:','LineWidth',2,...
                       'MarkerEdgeColor','r',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',10)
hold on
semilogy(4,relative_DR_y(4,3),'b^-','LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','b',...
                       'MarkerSize',10)
hold on
semilogy(1,relative_DR_y(1,5),'ko-.','LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','k',...
                       'MarkerSize',10)
hold on
legend('Correlated 0.4m','Correlated 0.7m','Correlated 1m','i.i.d','Fresnel')
%%%%%%
semilogy(1:Count_DR(1)-1,relative_DR_y(1:Count_DR(1)-1,1),'r:','LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',10)
hold on
semilogy(1:Count_DR(2)-1,relative_DR_y(1:Count_DR(2)-1,2),'r:','LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',10)
hold on
semilogy(1:Count_DR(4)-1,relative_DR_y(1:Count_DR(4)-1,4),'r:','LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',10)
hold on
semilogy(1:Count_DR(3)-1,relative_DR_y(1:Count_DR(3)-1,3),'b-','LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','b',...
                       'MarkerSize',10)
hold on
semilogy(1:150,relative_DR_y(1:150,5),'k-.','LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','k',...
                       'MarkerSize',10)
hold on

%%%%%%%%%


semilogy(1:5:Count_DR(1)-1,relative_DR_y(1:5:Count_DR(1)-1,1),'r^','LineWidth',2,...
                       'MarkerEdgeColor','r',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',10)
hold on
semilogy(2:5:Count_DR(2)-1,relative_DR_y(2:5:Count_DR(2)-1,2),'ro','LineWidth',2,...
                       'MarkerEdgeColor','r',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',10)
hold on
semilogy(3:5:Count_DR(4)-1,relative_DR_y(3:5:Count_DR(4)-1,4),'r*','LineWidth',2,...
                       'MarkerEdgeColor','r',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',10)
hold on
semilogy(4:5:Count_DR(3)-1,relative_DR_y(4:5:Count_DR(3)-1,3),'b^','LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','b',...
                       'MarkerSize',10)
hold on
semilogy(1:5:150,relative_DR_y(1:5:150,5),'ko','LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','k',...
                       'MarkerSize',10)
hold on





xlabel('Epoch')
ylabel('RE')

name=[Bigname,'/Error_pois_compare_5diff_type'];
h=gca;
set(h,'FontSize',14);
saveas(gcf,name,'png')
close

rate=zeros(5,1);
for i = 1:5
rate(i)=get_cr(log(relative_DR_y(1:Count_DR(i)-1,i)),15);



end
name_txt=[Bigname,'/convR.txt'];
fid=fopen(name_txt,'w');
fprintf(fid,'0.4m=%4.3s\n 0.7m=%4.3s\n iid=%4.3s\n 1m=%4.3s\n fresnel=%4.3s\n',...
    [rate(1) rate(2) rate(3) rate(4) rate(5)]);

fclose(fid);   


savepath=[Bigname,'/data.mat'];
save(savepath, 'relative_DR_y', 'rate');
