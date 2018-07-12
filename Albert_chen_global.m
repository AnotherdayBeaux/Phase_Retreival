%%% test albert chen global convergence

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
MaxIter=1000;
first_iter=2000;
Tol=1e-14;

SUBIM=256;

residual_DR_y=zeros(1,MaxIter);
relative_DR_x=zeros(1,MaxIter);
relative_DR_phase= zeros(1,MaxIter);


relative_DR2_x=zeros(1,MaxIter);
relative_DR2_phase= zeros(1,MaxIter);



num_im=1;

[IM,Y,Z,SNR,mask]=proc_(str,os_rate,num_mask,add_phase,...
                      sector_alpha,sector_beta,sector_mode,coer,SUBIM);
 if Noiseswitch == 0
     Z=abs(Y).^2;
 end
 B=Z;
 
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

lambda_0=rand(size(Y)).*exp(1i*rand(size(Y)))*2;

%%%
%%% add gaussian random noise
%NSR=0.05;
%a_r=NSR*norm(Y)/sqrt(size(Y,1)*size(Y,2));

%%%% end verify

    lambda_t=lambda_0;
    count_ini=1;  
    
    resi=1;
    AAlambda_pre=zeros(size(lambda_t));
    while resi > Tol && count_ini < MaxIter

            %%%%%% y_t+1 = argmin 1/2||y-A'A\lambda_t||^2+1/2gamma ||y-lambda_t||^2 %%%%%%%%
            Lam=zeros(size(B));
            abslam=abs(lambda_t);
            Lam(abslam==0)=1;
            lambda_tt=sqrt(B).*lambda_t./abs(lambda_t);
            lambda_tt(isnan(lambda_tt))=0;
            lambda_tt=lambda_tt+Lam.*sqrt(B);
            x_t=Nos_ifft(2*lambda_tt-lambda_t,os_rate,num_mask,mask);
            AAlambda=Nos_fft(x_t,os_rate,num_mask,mask);
            lambda_t=lambda_t+AAlambda-lambda_tt;
        
            ee=norm(IM(:)'*x_t(:),'fro')/(IM(:)'*x_t(:));
            rel_x=norm(ee*x_t-IM,'fro')/norm(IM,'fro');
            
            
            %ee=norm(Y(:)'*lambda_t(:),'fro')/(Y(:)'*lambda_t(:));
            ee=norm((Y(:)./abs(Y(:)))'*(lambda_t(:)./abs(lambda_t(:))),'fro')/((Y(:)./abs(Y(:)))'*(lambda_t(:)./abs(lambda_t(:))));
            rel_phase=norm(ee*(lambda_t./abs(lambda_t))-Y./abs(Y),'fro')/size(Y,1)/size(Y,2);
            
            
            
            resi=norm(abs(AAlambda)-abs(Y),'fro')/norm(abs(Y),'fro');
            
            residual_DR_y(count_ini)=resi;
            relative_DR_x(count_ini)=rel_x;
            relative_DR_phase(count_ini)=rel_phase;
            fprintf('count=%d\n resi=%f\n rel_x=%f\n rel_ph=%f\n', count_ini,resi,rel_x,rel_phase)
            count_ini=count_ini+1;
    
    end
    

phase_diff=relative_DR_phase(1:count_ini-1)-relative_DR_phase(2:count_ini);
find(phase_diff<0)

semilogy(1:count_ini-1, relative_DR_phase(1:count_ini-1),'r')



relx_diff=relative_DR_x(1:count_ini-1)-relative_DR_x(2:count_ini);
find(relx_diff<0)
hold on
semilogy(1:count_ini-1, relative_DR_x(1:count_ini-1),'g')



resi_diff=residual_DR_y(1:count_ini-1)-residual_DR_y(2:count_ini);
hold on
semilogy(1:count_ini-1, residual_DR_y(1:count_ini-1),'b')
find(resi_diff<0)    
legend('phase','rel_x','resi_y')
    
    
    
    
    
lambda_t=lambda_0;    
%%%%% DR on likelihood function.    
    count_DR=1;
    rel_x=1;
    while rel_x > Tol && count_DR< MaxIter
    
            x_t=Nos_ifft(lambda_t,os_rate,num_mask,mask); 
            
            AAlambda=Nos_fft(x_t,os_rate,num_mask,mask);
            y_tplus1=AAlambda;

            ylambda_k=2*y_tplus1-lambda_t;
            Q=abs(ylambda_k);
            % POISSON
            rot_z=(Q/gamma+sqrt(Q.^2/gamma^2+8*(2+1/gamma)*(Z)))/(4+2/gamma);
            
            % GAUSSIAN
            %rot_z=(sqrt(Z)+1/gamma*Q)/(1+1/gamma);
            
            z_tplus1=rot_z.*ylambda_k./Q;
            z_tplus1(isnan(z_tplus1))=0;
            QQ=zeros(size(Q));
            QQ(Q==0)=1;
            z_tplus1=z_tplus1+QQ.*sqrt(B);
    
            lambda_t=lambda_t+z_tplus1-y_tplus1;
   
            ee=norm(IM(:)'*x_t(:))/(IM(:)'*x_t(:));
            rel_x=norm(ee*x_t-IM,'fro')/norm(IM,'fro');
            relative_DR2_x(1,count_DR)=rel_x;
            
            ee=norm(Y(:)'*y_tplus1(:),'fro')/(Y(:)'*y_tplus1(:));
            rel_phase=norm(ee*y_tplus1./abs(y_tplus1)-Y./abs(Y),'fro')/size(Y,1)/size(Y,2);
            
            relative_DR2_phase(1,count_DR)=rel_phase;
            fprintf('count_DR=%d\n rel_x=%f\n rel_ph=%f\n', count_DR,rel_x,rel_phase);
            count_DR=count_DR+1;

     end
  
phase_diff=relative_DR2_phase(1:count_DR-1)-relative_DR2_phase(2:count_DR);
find(phase_diff<0)

semilogy(1:count_DR-1, relative_DR2_x(1:count_DR-1))
hold on 
semilogy(1:count_DR-1, relative_DR2_phase(1:count_DR-1),'r')
    