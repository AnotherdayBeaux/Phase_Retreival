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
            
print_im=0; 
ADMM=0;   % 0 admm not run; 1 admm run
% number of masks           
num_of_masks=floor(num_mask/2)+1;
                         
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

SUBIM=156;
[IM,Y,Z,SNR,mask]=proc_(str,os_rate,num_mask,add_phase,...
                      sector_alpha,sector_beta,sector_mode,coer,SUBIM);
 if Noiseswitch == 0
     Z=abs(Y).^2;
 end
 B= abs(Y);
 
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

% ADMM


mag_z_0=abs(rand(size(Y)));
phase_z_0=exp(0*2i*pi*rand(size(Y)));
z_0=mag_z_0.*phase_z_0;


lambda_0= (mag_z_0-B);
% DR



%3
%lambda_0=Y+rand(size(Y))*1e-1;




MaxIter=50000;
Tol=1e-10;
Gamma=0.8                                     %%%% to be modified

trial=0;



first_iter=1;

for gamma=Gamma
    relative_DR_y=zeros(length(Gamma),MaxIter);
    relative_ADMM_y=zeros(length(Gamma),MaxIter);
    %%%%%%%%%%%%%%
    trial=trial+1;
    %%%%%%%%%%%%%
    
    
    % ADMM
   
    
    lambda_t=lambda_0;
    z_t=z_0;
    
    
    count_ADMM=1;
    resi=1;
    %%%%%%%  initialization 

    
    z_0_DR=z_t;
    %%%%%%%
    while resi > Tol && count_ADMM < MaxIter
    
            %%%%%% y_t+1 = argmin 1/2||y-A'A\lambda_t||^2+1/2gamma ||y-lambda_t||^2 %%%%%%%%
            x_t=Nos_ifft(z_t-gamma*lambda_t,os_rate,num_mask,mask);
            AAlambda=Nos_fft(x_t,os_rate,num_mask,mask);
            %Resi=norm(y_tplus1,'fro')^2+norm(AAlambda,'fro')^2-2*abs(sum(sum(y_tplus1.*conj(AAlambda))));
            %resi1(count)=sqrt(norm(Resi))/norm(y_tplus1,'fro');
            y_tplus1=AAlambda;
    
            ylambda_k=y_tplus1+gamma*lambda_t;
            Q=abs(ylambda_k);
            % POISSON
            %rot_z=(Q/gamma+sqrt(Q.^2/gamma^2+8*(2+1/gamma)*B))/(4+2/gamma);
            
            % GAUSSIAN
            rot_z=(B+1/gamma*Q)/(1+1/gamma);
            ang_z_real=real(ylambda_k)./Q;
            ang_z_imag=imag(ylambda_k)./Q;

            z_t=rot_z.*ang_z_real+1i*rot_z.*ang_z_imag;
            %%%% lambda_t+1=lambda_t+z_t+1-y_t+1;
            
            lambda_t=lambda_t+(y_tplus1-z_t)/gamma;
            resi=norm(abs(z_t)-sqrt(Z),'fro')/norm(sqrt(Z),'fro');
            
            q_t=(abs(z_t)-gamma/(gamma+1)*B)*(gamma+1).*z_t./abs(z_t);
            Aq_t=Nos_ifft(q_t,os_rate,num_mask,mask);
            AAq_t=Nos_fft(Aq_t,os_rate,num_mask,mask);
            lambda_t_DR=2*AAq_t-q_t;
            
            
            ee=norm(Y(:)'*lambda_t_DR(:))/(Y(:)'*lambda_t_DR(:));
            rel_y=norm(ee*lambda_t_DR-Y,'fro')/norm(Y,'fro');
            relative_ADMM_y(trial, count_ADMM)=rel_y;
            fprintf('count_ADMM=%d\n resi=%f\n rel_y=%f\n', count_ADMM,resi, rel_y)
            count_ADMM=count_ADMM+1;
    end
    
    
    
    
    
    %%%%% DR
    mag_z_0=abs(z_0_DR);
    phase_z_0=mag_z_0./abs(z_0_DR);
    q_0_DR= (mag_z_0-gamma/(gamma+1)*B)*(gamma+1).*phase_z_0;
    Aq_0_DR=Nos_ifft(q_0_DR,os_rate,num_mask,mask);
    AAq_0_DR=Nos_fft(Aq_0_DR,os_rate,num_mask,mask);
    lambda_0_DR=2*AAq_0_DR-q_0_DR;
    lambda_t=lambda_0;
    count_DR=1;  
    
    resi=1;
 %gamma=1.2;
    while resi > Tol && count_DR < MaxIter
    
            %%%%%% y_t+1 = argmin 1/2||y-A'A\lambda_t||^2+1/2gamma ||y-lambda_t||^2 %%%%%%%%
            x_t=Nos_ifft(lambda_t,os_rate,num_mask,mask);
            AAlambda=Nos_fft(x_t,os_rate,num_mask,mask);
            %Resi=norm(y_tplus1,'fro')^2+norm(AAlambda,'fro')^2-2*abs(sum(sum(y_tplus1.*conj(AAlambda))));
            %resi1(count)=sqrt(norm(Resi))/norm(y_tplus1,'fro');
            y_tplus1=AAlambda;
    
            ylambda_k=2*y_tplus1-lambda_t;
            Q=abs(ylambda_k);
            % POISSON
            %rot_z=(Q/gamma+sqrt(Q.^2/gamma^2+8*(2+1/gamma)*B))/(4+2/gamma);
            
            % GAUSSIAN
            rot_z=(B+1/gamma*Q)/(1+1/gamma);
            ang_z_real=real(ylambda_k)./Q;
            ang_z_imag=imag(ylambda_k)./Q;

            z_tplus1=rot_z.*ang_z_real+1i*rot_z.*ang_z_imag;
    
            %%%% lambda_t+1=lambda_t+z_t+1-y_t+1;
            lambda_t=lambda_t+z_tplus1-y_tplus1;
            resi=norm(abs(y_tplus1)-sqrt(Z),'fro')/norm(sqrt(Z),'fro');
            ee=norm(Y(:)'*lambda_t(:))/(Y(:)'*lambda_t(:));
            rel_y=norm(ee*lambda_t-Y,'fro')/norm(Y,'fro');
            relative_DR_y(trial, count_DR)=rel_y;
            fprintf('count_DR=%d\n resi=%f\n rel_y=%f\n', count_DR,resi, rel_y)
            %fprintf('count=%d\n lambda_2=%f\n resi=%f\n', count,lambda_2,resi)
            count_DR=count_DR+1;
    
    end
    
    
    %%%%% semilogy admm/dr
    figure
    semilogy(1:count_DR-1,relative_DR_y(trial,1:count_DR-1),'g')
    hold on
    semilogy(1:count_ADMM-1,relative_ADMM_y(trial,1:count_ADMM-1),'r')
    title(['gamma= ',num2str(gamma)])
    xlabel('Iter')
    ylabel('relative error')
    legend('DR','ADMM')
    %%%%%%%%%%%%%%%%%%%%%%%

end


%%%%% compare different gamma in the DR/ADMM method