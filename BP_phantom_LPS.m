%%%% Blind  ptychography  no-noise case

%%%% |phi(n)-phi_0| < \delta * \pi; phi_0 is uniform for all element in the
%%%% mask

%%%% The main task of this program is to
%%%%     1. test  Phantom case  (find optimal linear phase)
%%%%     2. add salt noise
%%%%     3. enlarge the size of mask to incorporate the support of image.
%%%%      modify on 6/6/2018




Bigname='/Users/Beaux/Desktop/phase retrieval'   ; % location of image

str1=[Bigname,'/image_lib/cameraman.bmp'];   %% consider Im= Cameraman+i* barbera
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

         
MaxIter=300;
Tol=1e-13;

SUBIM=256;   % truncate image size cameraman default 256-by-256

bd_con=1;     % periodic boundary condition.
p_c=0;        % perturb center

type_im=1;    % type_im=0: real im str1
              % type_im=1; str1+i*str1
              % type_im=2; str1+i*str2
              % type_im=3; str1* rand_phase
mask_type=0;  % 2 fresnel mask
              % 1 correlated mask
              % 0 iid mask
coer=1;

salt_on=1;
salt_ratio=0.00;
gamma=1;
relative_DR_y=zeros(MaxIter,1);
relative_DR_maskLPS =zeros(MaxIter,1);
relative_DR_yshift=zeros(MaxIter,5);
resi_DR_y=zeros(MaxIter,1);
o_l = 0.5;%0.2:0.1:0.8 ; % overlaping
cim_diff=30; %20:15:120; % make sure each mask shot overlaps support of phantom
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
    Bigname='/Users/Beaux/Desktop/phase retrieval/code/16/LPS_CamBar';
    mkdir(Bigname)
    
   
    
    [IM,salt_noise, Y,Z,NSR,phase_arg,nor_ptycho,nor_phase]=ini_rand_ptycho_salt(str1,str2,os_rate,num_mask,coer,SUBIM,...
    type_im,l_patch,x_c_p,y_c_p,bd_con,mask_type,c_l, beta_1,beta_2, rho,salt_on,salt_ratio);
 
 
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
    
    
    
    
    for tr= 1:5
    
    
    
    
    k1=-1;%floor(rand(1)*n/m/4);
    l1=1;%floor(rand(1)*n/m/4);
    a1=1*exp(2i*pi*rand(1));
    
    k2=1;
    l2=0;a2=1*exp(2i*pi*rand(1));
    
    k3=0;
    l3=1;a3=1*exp(2i*pi*rand(1));
    
    k4=1;
    l4=1;a4=1*exp(2i*pi*rand(1));
    
    k5=1;
    l5=2;a5=1*exp(2i*pi*rand(1));
    
    
    
    im_x_grid=(1:SUBIM)'*ones(1,SUBIM);
    im_y_grid=ones(SUBIM,1)*(1:SUBIM);
    %l_phase_drift_im=a*exp(2i*pi/n*(im_x_grid*k+im_y_grid*l))+...
    %                a1*exp(2i*pi/n*(im_x_grid*k1+im_y_grid*l1))+...
    %                a2*exp(2i*pi/n*(im_x_grid*k2+im_y_grid*l2))+...
    %                a3*exp(2i*pi/n*(im_x_grid*k3+im_y_grid*l3))+...
    %                a4*exp(2i*pi/n*(im_x_grid*k4+im_y_grid*l4))   ;
    %l_phase_drift_im = l_phase_drift_im ./ abs(l_phase_drift_im) ;     
    l_phase_drift_im  = @(kk,ll) exp(2i*pi/n*(im_x_grid*kk+im_y_grid*ll));
    
    
    mask_x_grid=(1:l_patch)'*ones(1,l_patch);
    mask_y_grid=ones(l_patch,1)*(1:l_patch);
    l_phase_drift_mask= @(kk,ll) exp(-2i*pi/n*(mask_x_grid*kk+mask_y_grid*ll));
    
    %%% IN THIS PROGRAM, we apply sum (a_j e^{ikr})
    l_phase_drift_mask_sum = a1*l_phase_drift_mask(k1,l1)+...
                         a2*l_phase_drift_mask(k2,l2)+...
                         a3*l_phase_drift_mask(k3,l3)+...
                         a4*l_phase_drift_mask(k4,l4)+...
                         a5*l_phase_drift_mask(k5,l5)    ;
                     
    l_phase_drift_mask_sum = l_phase_drift_mask_sum./abs(l_phase_drift_mask_sum);                 
    %%%
    
    
    lambda_0=10*(rand(size(Y))+1i*rand(size(Y))); % random initialize  image
    
    est_phase=delta*(rand(l_patch,l_patch)-1/2);
    mask=exp(2i*pi*(phase_arg));
    %DR_Pois
    
    
    lambda_t=lambda_0;
    
    %%%%% correlated perturbation.
    
    
    %mask_estimate= exp(2i*pi*phase_arg).*l_phase_drift_mask.*  exp(1/2*2i*pi*(rand(l_patch,l_patch)-1/2));
    %mask_estimate= exp(2i*pi*phase_arg).*l_phase_drift_mask.*  Cor_Mask(m,c_l,1/2);
    
    mask_estimate= exp(2i*pi*phase_arg).*l_phase_drift_mask_sum.* exp(1/11*2i*pi*(rand(l_patch,l_patch)-1/2)) ;
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
    res_x=1;

    
    
    
    
    
    mask_estimate_ini=mask_estimate;
while res_x > Tol && count_DR< MaxIter
     MaxIter_u_ip = 30;
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
            % POISSON
            rot_z=(Q/gamma+sqrt(Q.^2/gamma^2+8*(2+1/gamma)*(Z)))/(4+2/gamma);
            
            %GAUSSIAN
            %rot_z=(sqrt(Z)+1/gamma*Q)/(1+1/gamma);
            
            z_tplus1=rot_z.*ylambda_k./Q;
            % test if lambda_t converge to a fixed point
            lambda_tpre=lambda_t;
           
            lambda_t=lambda_t+z_tplus1-y_tplus1;
            rel_q=norm(lambda_t-lambda_tpre,'fro')/norm(lambda_t,'fro');
           
            % end test
            ee=norm(IM(:)'*x__t(:))/(IM(:)'*x__t(:));
            
            rel_xx=norm(ee*x__t-IM,'fro')/norm(IM,'fro');
            res_x= norm(abs(Y)-abs(y_tplus1),'fro')/norm_Y;
            %fprintf('u_i=%d\n rel_x=%f\n res_x=%f\n', update_im,rel_xx,res_x);
            %fprintf('u_i=%d\n rel_q=%f rel_x=%f res_x=%f\n', update_im, rel_q,rel_x,res_x);
            residual_x(update_im+3)=res_x;
            resi_diff_im = 1/3*norm(residual_x(update_im+1:update_im+3)-residual_x(update_im:update_im+2),1)/res_x;
            update_im=update_im+1;
            x__tpre=x__t;
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
            %%%% POISSON
            rot_z=(Q/gamma+sqrt(Q.^2/gamma^2+8*(2+1/gamma)*(Z)))/(4+2/gamma);
            
            %%%% GAUSSIAN
            %rot_z=(sqrt(Z)+1/gamma*Q)/(1+1/gamma);
            
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
    
    
    
    mask_t=mask_estimate./abs(mask_estimate);
    %%%% MPC enforcement.
    %mask_arg=angle(mask_estimate./mask);       % center
    %center_arg=sum(sum(mask_arg))/m/m;
    %proj_mask_arg=max(mask_arg-center_arg, -delta*pi);
    %proj_mask_arg=min(proj_mask_arg, delta*pi);
    %mask_estimate = exp(1i*proj_mask_arg);
    %%%% MPC enforcement.
    
    %shift_IM_1=IM.*l_phase_drift_im(k1,l1);
   
    
     
    %%%%% min || IM- exp(ikn) x_t||
    ee=norm(IM(:)'*x__t(:))/(IM(:)'*x__t(:));
    [rec_k1,rec_l1]= getrid_LPS_m(IM, ee*x__t,n);
    LPS_x_t=exp(-2i*pi/n*(im_x_grid*rec_k1+im_y_grid*rec_l1)).*x__t;
    ee=norm(IM(:)'*LPS_x_t(:))/(IM(:)'*LPS_x_t(:));
    rel_LPS_x=norm(LPS_x_t*ee-IM,'fro')/norm(IM,'fro');
    
     %%%%% min || mask- exp(ikm) mask_estimate||
    ee_mask=norm(mask(:)'*mask_t(:))/(mask(:)'*mask_t(:));
    [recm_k1,recm_l1]= getrid_LPS_m(mask, ee_mask*mask_t,n);
    mask_t_LPS= exp(-2i*pi/n*(mask_x_grid*recm_k1+mask_y_grid*recm_l1)).*mask_t;
    ee_mask=norm(mask(:)'*mask_t_LPS(:))/(mask(:)'*mask_t_LPS(:));
    rel_LPS_mask= norm(ee_mask*mask_t_LPS-mask,'fro')/norm(mask,'fro');
    
    
    resi_DR_y(count_DR,1)= res_x;
    relative_DR_y(count_DR,1)=rel_LPS_x;
    relative_DR_maskLPS(count_DR,1)=rel_LPS_mask;
    
    fprintf('count_DR=%d\n rel_LPS_x=%4.3s rel_LPS_mask=%4.3s res_x=%4.3s\n',...
                count_DR,rel_LPS_x, rel_LPS_mask, res_x);
    count_DR=count_DR+1;
    
    
   
end 



%%%%% save the convergence rate



last=15;
rate_resi=get_cr(log(resi_DR_y(1:count_DR-1,1)),last);
rate_rel_im=get_cr(log(relative_DR_y(1:count_DR-1)),last);
rate_rel_mask=get_cr(log(relative_DR_maskLPS(1:count_DR-1)),last);

foldername=[Bigname,'/trial=',num2str(tr),', mask size = ',num2str(l_patch)];
mkdir(foldername)
name_txt=[foldername,'/convR.txt'];
fid=fopen(name_txt,'w');
fprintf(fid,'resiR=%s  relRim=%s relRmask=%s\n imk1=%f iml1=%f \n maskk1=%f maskl1=%f\n',...
    [rate_resi rate_rel_im rate_rel_mask rec_k1 rec_l1 recm_k1 recm_l1]);

fclose(fid);
%%%%%%%% end save

figure
semilogy(1,resi_DR_y(1,1),'r*:','LineWidth',2,...
                       'MarkerEdgeColor','r',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',12)
hold on
semilogy(2,relative_DR_y(2,1),'bo-.','LineWidth',2,...
                       'MarkerEdgeColor','b',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',12)
hold on
semilogy(3,relative_DR_maskLPS(3,1),'k^-.','LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',12)

xlabel('Epoch')
ylabel(' ')
legend('Residual','RE image','RE mask')


semilogy(1:count_DR-1,resi_DR_y(1:count_DR-1,1),'r:','LineWidth',2,...
                       'MarkerEdgeColor','r',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',12)
hold on
semilogy(1:count_DR-1,relative_DR_y(1:count_DR-1,1),'b-.','LineWidth',2,...
                       'MarkerEdgeColor','b',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',12)
hold on
semilogy(1:count_DR-1,relative_DR_maskLPS(1:count_DR-1,1),'k-.','LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',12)
hold on

%%%%%%%%
semilogy(1:10:count_DR-1,resi_DR_y(1:10:count_DR-1,1),'r*','LineWidth',2,...
                       'MarkerEdgeColor','r',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',12)
hold on
semilogy(2:10:count_DR-1,relative_DR_y(2:10:count_DR-1,1),'bo','LineWidth',2,...
                       'MarkerEdgeColor','b',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',12)
hold on
semilogy(3:10:count_DR-1,relative_DR_maskLPS(3:10:count_DR-1,1),'k^','LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',12)


foldername=[Bigname,'/trial=',num2str(tr),', mask size = ',num2str(l_patch)];
name=[foldername,'/convergence_LPS_phantom'];
grid on
h=gca;
set(h,'FontSize',14);
saveas(gcf,name,'png')
close
savepath=[foldername,'/data.mat'];
save(savepath, 'relative_DR_y', 'resi_DR_y','relative_DR_maskLPS','count_DR');

%%%%%%%%%
figure
imshow(abs(LPS_x_t*ee)/255);
title('magnitude reconstruction')
h=gca;
set(h,'FontSize',14);
name=[foldername,'/modular_phantom_recon'];
saveas(gcf,name,'png')
close


figure
imshow(abs(IM)/255);
title('magnitude phantom image')
h=gca;
set(h,'FontSize',14);
name=[foldername,'/modular_phantom_ini_salt'];
saveas(gcf,name,'png')
close


figure
ang_phantom=angle(IM);
colormap('hot')
imagesc(ang_phantom);
colorbar
title('phantom phase')
h=gca;
set(h,'FontSize',14);
name=[foldername,'/phase_phantom_salt'];
saveas(gcf,name,'png')
close
%figure
%shift_mask=ee_mask*mask_t_LPS;
%diff_angle_mask=angle(mask.*conj(shift_mask));
%colormap('hot')
%imagesc(diff_angle_mask)
%colorbar
%name=[foldername,'/mask_recon_error'];
%saveas(gcf,name,'png')
%close





figure
LPS_x_t_m=exp(2i*pi/n*(im_x_grid*recm_k1+im_y_grid*recm_l1)).*x__t;
ee_m=norm(IM(:)'*LPS_x_t_m(:))/(IM(:)'*LPS_x_t_m(:));
ang_im_diff=angle(ee_m*LPS_x_t_m.*conj(IM));
colormap('hot')
imagesc(ang_im_diff)
colorbar
name=[foldername,'/im_phase_error'];
h=gca;
set(h,'FontSize',14);
saveas(gcf,name,'png')
close


figure
ang_mask_diff=angle(ee_mask*mask_t_LPS.*conj(mask));
colormap('hot')
imagesc(ang_mask_diff)
colorbar
name=[foldername,'/mask_phase_error'];
h=gca;
set(h,'FontSize',14);
saveas(gcf,name,'png')
close








    
%save_ws=[foldername,'/all_variable.mat'];
%save(save_ws)

    end
    
    
    
   