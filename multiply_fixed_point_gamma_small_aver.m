%%%% support my idea that when rho is small the basin of attraction might
%%%% be multiple. 



Bigname='/Users/Beaux/Desktop/phase retrieval'   ; % location of image


str=[Bigname,'/image_lib/Cameraman.bmp'];

os_rate=2;  %% oversampling ratio
TVswitch=0;    % TV   1==ON
            %      0==OFF
tv=0.1; % tv * ||grad_image||_1
Noiseswitch = 0; % 1= add poisson noise to |Y|^2;        
num_mask=3; % 0  0 mask / Id mask
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
            
coer=1;    % regulating the SNR. coer higher, noise higher.
MaxIter=20000;
Tol=1e-14;
SUB_IM=ones(1,1)*256;
gamma=0.2;
l1=length(SUB_IM);

residual_DR_y=zeros(l1,MaxIter);
residual_DR_y_noini=zeros(l1,MaxIter);
residual_DR_y_noini_relative=zeros(l1,MaxIter);
num_im=1;
for SUBIM=SUB_IM

[IM,Y,Z,SNR,mask]=proc_(str,os_rate,num_mask,add_phase,...
                      sector_alpha,sector_beta,sector_mode,coer,SUBIM);
 if Noiseswitch == 0
     Z=abs(Y).^2;
 end
 
 
[Na,Nb]=size(IM);
norm_IM=norm(IM,'fro');
n=Na*Nb; %% dim of data vector
[Y_Na,Y_Nb]=size(Y); %%% This Y_Nb= Nb*os_rate*'num_mask'
% Q=ifft2(Y(1:os_rate:os_rate*Na,1:os_rate:os_rate*Nb));
% a=norm((Q-IM),'fro');

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

% DR
lambda_0=rand(size(Y)).*Y.*exp(1i*1e-2*rand(size(Y)))*2;
%lambda_0=abs(Y);
% ADMM
%3
%lambda_0=Y+rand(size(Y))*1e-1;
lambda_2=0;




    lambda_t=lambda_0;
    count_ini=1;  
    
    resi=1;
    while resi > Tol && count_ini < MaxIter
           
            %%%%%% y_t+1 = argmin 1/2||y-A'A\lambda_t||^2+1/2gamma ||y-lambda_t||^2 %%%%%%%%
            x_t=Nos_ifft(lambda_t,os_rate,num_mask,mask);
            AAlambda=Nos_fft(x_t,os_rate,num_mask,mask);
            %Resi=norm(y_tplus1,'fro')^2+norm(AAlambda,'fro')^2-2*abs(sum(sum(y_tplus1.*conj(AAlambda))));
            %resi1(count)=sqrt(norm(Resi))/norm(y_tplus1,'fro');
            y_tplus1=AAlambda;
    
            %%%%% z_t+1 = argmin |z|^2-b log(|z|^2)+1/2gamma*||2y_t+1-lambda_t-z||^2
            ylambda_k=2*y_tplus1-lambda_t;
            Q=abs(ylambda_k);
            %Poisson likelihood 
            rot_z=(Q/gamma+ sqrt(Q.^2/gamma^2+8*(2+1/gamma)*Z))/(4+2/gamma);
            %Gaussian likelihood
            %rot_z=(sqrt(Z)+1/gamma*Q)/(1+1/gamma);
            ang_z_real=real(ylambda_k)./Q;
            ang_z_imag=imag(ylambda_k)./Q;

            z_tplus1=rot_z.*ang_z_real+1i*rot_z.*ang_z_imag;
    
            %%%% lambda_t+1=lambda_t+z_t+1-y_t+1;
            %lambda_2=0*norm(y_tplus1-lambda_t)/norm(y_tplus1);
            
            %Lambda_2(count)=lambda_2;
            %lambda_2=lambda_t-AAlambda;
            %resi1=norm(lambda_2,'fro')/norm(lambda_t,'fro');
            lambda_t=lambda_t+z_tplus1-y_tplus1;
    
            resi=norm(abs(y_tplus1)-sqrt(Z),'fro')/norm(sqrt(Z),'fro');
            residual_DR_y(num_im, count_ini)=resi;
            
            
            
            %residual_DR_lambda(num_im, count)=resi1;
            %fprintf('count=%d\n norm(lambda_2)=%f\n resi=%f\n',count,resi1,resi)
            
            %IM_resi=norm(x_t,'fro')^2+norm_IM^2-2* abs(sum(sum(IM.*conj(x_t))));
            %dev(count)=sqrt(norm(IM_resi))/norm_IM;
            fprintf('count_ini=%d\n resi=%s\n', count_ini,resi)
            count_ini=count_ini+1;
            %figure
            %imshow(abs(x_t)/255)
            %pause(0.5)
    
    end
    
    
    
    
    
    
    
    
    %%%%%%%  no initialization 
    lambda_0=rand(size(Y)).*Y.*exp(1i*2*pi*rand(size(Y)));
    lambda_t=lambda_0;
    y_t=lambda_t;
    count_noini=1;  
    
    resi=1;
 %for ii =1:3
    while resi > Tol && count_noini < MaxIter
            
            %%%%%% y_t+1 = argmin 1/2||y-A'A\lambda_t||^2+1/2gamma ||y-lambda_t||^2 %%%%%%%%
            x_t=Nos_ifft(lambda_t,os_rate,num_mask,mask);
            AAlambda=Nos_fft(x_t,os_rate,num_mask,mask);
            %Resi=norm(y_tplus1,'fro')^2+norm(AAlambda,'fro')^2-2*abs(sum(sum(y_tplus1.*conj(AAlambda))));
            %resi1(count)=sqrt(norm(Resi))/norm(y_tplus1,'fro');
            y_tplus1=AAlambda;
    
            %%%%% z_t+1 = argmin |z|^2-b log(|z|^2)+1/2gamma*||2y_t+1-lambda_t-z||^2
            ylambda_k=2*y_tplus1-lambda_t;
            Q=abs(ylambda_k);
            %Poisson likelihood 
            %rot_z=(Q/gamma+ sqrt(Q.^2/gamma^2+8*(2+1/gamma)*Z))/(4+2/gamma);
            %Gaussian likelihood
            rot_z=(sqrt(Z)+1/gamma*Q)/(1+1/gamma);
            ang_z_real=real(ylambda_k)./Q;
            ang_z_imag=imag(ylambda_k)./Q;

            z_tplus1=rot_z.*ang_z_real+1i*rot_z.*ang_z_imag;
    
            %%%% lambda_t+1=lambda_t+z_t+1-y_t+1;
            %lambda_2=0*norm(y_tplus1-lambda_t)/norm(y_tplus1);
            
            %Lambda_2(count)=lambda_2;
            %lambda_2=lambda_t-AAlambda;
            %resi1=norm(lambda_2,'fro')/norm(lambda_t,'fro');
            lambda_t=lambda_t+z_tplus1-y_tplus1;
            
            resi=norm(abs(y_tplus1)-sqrt(Z),'fro')/norm(sqrt(Z),'fro');
            residual_DR_y_noini(num_im, count_noini)=resi;
            
            %%% calculate relative error
            ee=norm(y_t(:)'*y_tplus1(:))/(y_t(:)'*y_tplus1(:));
            rel_y=norm(ee*y_tplus1-y_t,'fro')/norm(y_t,'fro');
            residual_DR_y_noini_relative(num_im,count_noini)=rel_y;
            y_t=y_tplus1;
            
            %residual_DR_lambda(num_im, count)=resi1;
            %fprintf('count=%d\n norm(lambda_2)=%f\n resi=%f\n',count,resi1,resi)
            
            %IM_resi=norm(x_t,'fro')^2+norm_IM^2-2* abs(sum(sum(IM.*conj(x_t))));
            %dev(count)=sqrt(norm(IM_resi))/norm_IM;
            fprintf('count_noini=%d\n resi=%s\n relative=%s\n', count_noini,resi,rel_y)
            count_noini=count_noini+1;
            %figure
            %imshow(abs(x_t)/255)
            %pause(0.5)
    
    end
 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
  num_im=num_im+1  
  
  
  
  
  
  
  
  
    
end



name='/Users/Beaux/Desktop/phase retrieval/code/report 3.8';
figure
semilogy(1:3000,mean(residual_DR_y(:,1:3000)),'g')

xlabel('iteration')
ylabel('residual')
title('semilogy of convergence, with knowledge of phase (averaged)')
Sname=[name,'/rho_large_knowledge_of_phase_p'];
saveas(gcf,Sname,'png')

figure
semilogy(1:20000,mean(residual_DR_y_noini(:,1:20000)),'g')
xlabel('iteration')
ylabel('residual')
title('semilogy of convergence, without knowledge of phase (averaged)')    

Sname=[name,'/rho_large_no_knowledge_of_phase_residual_p'];
saveas(gcf,Sname,'png')
    
figure
semilogy(1:20000,mean(residual_DR_y_noini_relative(:,1:20000)),'g')
xlabel('iteration')
ylabel('relative')
title('semilogy of convergence relative, without knowledge of phase (averaged)')    

Sname=[name,'/rho_large_no_knowledge_of_phase_relative_p'];
saveas(gcf,Sname,'png')

