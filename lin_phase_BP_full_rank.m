%%%% Blind  ptychography  no-noise case

%%%% |phi(n)-phi_0| < \delta * \pi; phi_0 is uniform for all element in the
%%%% mask

%%%% The main task of this program is to
%%%%     1. test the dependence of success
%%%%        rate on \delta  (full rank perturb, rand mask ini)





Bigname='/Users/Beaux/Desktop/phase retrieval'   ; % location of image

str1=[Bigname,'/image_lib/Cameraman.bmp'];   %% consider Im= Cameraman+i* barbera
str2=[Bigname,'/image_lib/barbara.png'];

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

         
MaxIter=100;
Tol=1e-8;

SUBIM=156;   % truncate image size cameraman default 256-by-256

bd_con=1;     % periodic boundary condition.
p_c=0;        % perturb center
type_im=2;    % str1+i*str2 
mask_type=2;  % 2 fresnel mask
              % 1 correlated mask
              % 0 iid mask
coer=1;



gamma=1;
relative_DR_y=zeros(MaxIter,1);
relative_DR_yshift=zeros(MaxIter,1);
resi_DR_y=zeros(MaxIter,1);
o_l = 0.5;%0.2:0.1:0.8 ; % overlaping
cim_diff=30; %20:15:120;
perturb = 3;    % set the max magnitude of either horizontal or vertical perturb
rank_pert = 100;% =1 rank_one perturb
  % delta < 1/2 regulates the MPC.
delta=1/2;
cor_m=1;

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
     c_l=l_patch;


    %Bigname = '/Users/Beaux/Desktop/phase retrieval/code/16/linear phase/full_rank_randphased_shift_0.5';
    %mkdir(Bigname);
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
    Bigname=['/Users/Beaux/Desktop/phase retrieval/code/16/Ini_FresnelMask_q=',num2str(q)];
    mkdir(Bigname)
    
    count_rho=1;
    for rho = 0.5
   
    
    [IM,Y,Z,NSR,phase_arg,nor_ptycho,nor_phase]=ini_rand_ptycho(str1,str2,os_rate,num_mask,coer,SUBIM,...
    type_im,l_patch,x_c_p,y_c_p,bd_con,mask_type,c_l, delta, beta_1,beta_2, rho);
 
 
    [Na,Nb]=size(IM);
    norm_IM=norm(IM,'fro');

    [Y_Na,Y_Nb]=size(Y);
    
    Z=abs(Y).^2;
    b=abs(Y);
    norm_Y=norm(b,'fro');
    B=sqrt(abs(Z));  %%% noise case
    
    %%% linear phase shift generator
    m=l_patch;
    n=SUBIM;
    k=-1;%floor(rand(1)*n/m/4);
    l=1;%floor(rand(1)*n/m/4);
    k1=0;
    l1=0;
    im_x_grid=(1:SUBIM)'*ones(1,SUBIM);
    im_y_grid=ones(SUBIM,1)*(1:SUBIM);
    l_phase_drift_im=exp(2i*pi/n*(im_x_grid*(k+k1)+im_y_grid*(l+l1)));
    
    mask_x_grid=(1:l_patch)'*ones(1,l_patch);
    mask_y_grid=ones(l_patch,1)*(1:l_patch);
    l_phase_drift_mask = exp(-2i*pi/n*(mask_x_grid*(k+k1)+mask_y_grid*(l+l1)));
    %%%
    
    
    lambda_0=10*(rand(size(Y))+1i*rand(size(Y))); % random initialize  image
    
    est_phase=delta*(rand(l_patch,l_patch)-1/2);
    mask=exp(2i*pi*(phase_arg));
    %DR_Pois
    
    
    lambda_t=lambda_0;
    
    %%%%% correlated perturbation.
    
    
    %mask_estimate= exp(2i*pi*phase_arg).*l_phase_drift_mask.*  exp(1/2*2i*pi*(rand(l_patch,l_patch)-1/2));
    %mask_estimate= exp(2i*pi*phase_arg).*l_phase_drift_mask.*  Cor_Mask(m,c_l,1/2);
    
    %mask_estimate= exp(2i*pi*phase_arg).*l_phase_drift_mask.* exp(1/2*2i*pi*(rand(l_patch,l_patch)-1/2)) ;
    
    phase_arg_estimate= 1/2 * 0.7/l_patch *( (((1:l_patch)-13).^2)'*ones(1,l_patch)+...
            (ones(l_patch,1)*((1:l_patch)-14).^2)   );  % argument of phase
    mask_estimate = exp(2i*pi*phase_arg_estimate).* exp(1/5*2i*pi*(rand(l_patch,l_patch)-1/2)) ;
    
    
    
    
    
    
    
    
    
    
    %%%%%
    
    
    
    %%%% MPC enforcement.
    %mask_arg=angle(mask_estimate./mask);       % center
    %center_arg=sum(sum(mask_arg))/m/m;
    %proj_mask_arg=max(mask_arg-center_arg, -delta*pi);
    %proj_mask_arg=min(proj_mask_arg, delta*pi);
    %mask_estimate = exp(1i*proj_mask_arg);
    %%%% MPC enforcement.
    
    
    
    
    
    %enforce_mask = mask_estimate./exp(2i*pi*phase_arg-1i*pi/2) ;
    %enf_imag=max(imag(enforce_mask),0)+1e-5;
    %mask_estimate = (real(enforce_mask)+1i*enf_imag).*exp(2i*pi*phase_arg-1i*pi/2);
    %mask_estimate= mask_estimate./abs(mask_estimate);
    
    %mask_estimate= exp(2i*pi*(est_phase));
    %lambda_t=ptycho_rand_fft(IM.*l_phase_drift_im,os_rate,num_mask,mask_estimate,l_patch,x_c_p,y_c_p,bd_con);
    %lambda_t=0+10*rand(size(lambda_t)).*exp(2*1i*pi*rand(size(lambda_t)));
    
    count_DR=1;   
    rel_x=1;
    shift_IM=IM.*l_phase_drift_im;
    mask_estimate_ini=mask_estimate;
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
    while  update_im < 5 || ( update_im < MaxIter_u_ip && resi_diff_im >1e-4)
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
            
            z_tplus1=rot_z.*ylambda_k./Q;
         
            lambda_t=lambda_t+z_tplus1-y_tplus1;
            ee=norm(IM(:)'*x__t(:))/(IM(:)'*x__t(:));

            rel_xx=norm(ee*x__t-IM,'fro')/norm(IM,'fro');
            res_x= norm(abs(Y)-abs(y_tplus1),'fro')/norm_Y;
            %fprintf('u_i=%d\n rel_x=%f\n res_x=%f\n', update_im,rel_xx,res_x);
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
            rot_z=(Q/gamma+sqrt(Q.^2/gamma^2+8*(2+1/gamma)*(Z)))/(4+2/gamma);
            z_mask_tplus1 = rot_z.* Q_phase./Q;
            
            
            fft_phase=fft_phase+ z_mask_tplus1   -fft_phase_1;
            
            
            ee=norm(mask(:)'*mask_tplus1(:))/(mask(:)'*mask_tplus1(:));
            rel_mask=norm(ee*mask_tplus1 - mask,'fro')/norm(mask,'fro');
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
    eee=norm(shift_IM(:)'*x__t(:))/(shift_IM(:)'*x__t(:));
    rel_shift_x=norm(eee*x__t-shift_IM,'fro')/norm(shift_IM,'fro');
    rel_x=norm(ee*x__t-IM,'fro')/norm(IM,'fro');
    
    resi_DR_y(1,count_DR)= res_x;
    relative_DR_y(1,count_DR)=rel_x;
    relative_DR_yshift(1,count_DR)=rel_shift_x;
    fprintf('count_DR=%d\n rel_x=%s res_x=%s\n', count_DR,rel_x, res_x);
    count_DR=count_DR+1;
    
    
   
end 
foldername=[Bigname,'/rho=',num2str(rho)];
mkdir(foldername)


figure
semilogy(1:count_DR-1,relative_DR_y(1,1:count_DR-1),'g')
hold on 
semilogy(1:count_DR-1,relative_DR_yshift(1,1:count_DR-1),'g.')
xlabel('iter')
ylabel(' ')
hold on
semilogy(1:count_DR-1,resi_DR_y(1,1:count_DR-1),'r')
legend('rel err x_0','rel err linear-shift x_0','residual')
name=[foldername,'/convergence'];
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
    
    
    
    end_error(count_rho)=rel_x;
    end_error_resi(count_rho)=res_x;
    count_rho=count_rho+1;
    end
    
    semilogy(0.2:0.2:2, end_error, '-bs')
    hold on 
    semilogy(0.2:0.2:2, end_error_resi, '-rs')
    xlabel('rho=')
    legend('Rel Err', 'Residual')
    
    
    figure
    colormap('hot')
    fft2_IM=fft2(IM);
    imagesc(abs(fft2_IM))
    colorbar
    title('fft of original IM')
    
    
    
    
    figure 
    colormap('hot')
    fft2_reconIM=fft2(x__t);
    imagesc(abs(fft2_reconIM))
    colorbar
    title('fft of recon IM')
    
    
    
    
    
    