Bigname='/Users/Beaux/Desktop/phase retrieval'   ; % location of image


str=[Bigname,'/image_lib/barbara.png'];

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

SUB_IM=156:20:512;
gamma=1.5;
l1=length(SUB_IM);

residual_DR_y=zeros(l1,MaxIter);
residual_DR_lambda=zeros(l1,MaxIter);
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
%lambda_0=rand(size(Y)).*Y.*exp(1i*1e-1*rand(size(Y)))*10;
lambda_0=abs(Y);
% ADMM
%3
%lambda_0=Y+rand(size(Y))*1e-1;
lambda_2=0;




MaxIter=1000;
Tol=1e-10;

    %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%
    lambda_t=lambda_0;
    count=1;  
    
    resi=1;
    resi1=inf;
    tic;
 %gamma=1.2;
 lambda_t_1=lambda_t;
 %for ii =1:3
    while resi > Tol && count < MaxIter
            close
            %%%%%% y_t+1 = argmin 1/2||y-A'A\lambda_t||^2+1/2gamma ||y-lambda_t||^2 %%%%%%%%
            x_t=Nos_ifft(lambda_t,os_rate,num_mask,mask);
            AAlambda=Nos_fft(x_t,os_rate,num_mask,mask);
            %Resi=norm(y_tplus1,'fro')^2+norm(AAlambda,'fro')^2-2*abs(sum(sum(y_tplus1.*conj(AAlambda))));
            %resi1(count)=sqrt(norm(Resi))/norm(y_tplus1,'fro');
            y_tplus1=AAlambda;
    
            %%%%% z_t+1 = argmin |z|^2-b log(|z|^2)+1/2gamma*||2y_t+1-lambda_t-z||^2
            ylambda_k=2*y_tplus1-lambda_t;
            Q=abs(ylambda_k);
            rot_z=(Q/gamma+ sqrt(Q.^2/gamma^2+8*(2+1/gamma)*Z))/(4+2/gamma);
            ang_z_real=real(ylambda_k)./Q;
            ang_z_imag=imag(ylambda_k)./Q;

            z_tplus1=rot_z.*ang_z_real+1i*rot_z.*ang_z_imag;
    
            %%%% lambda_t+1=lambda_t+z_t+1-y_t+1;
            %lambda_2=0*norm(y_tplus1-lambda_t)/norm(y_tplus1);
            
            %Lambda_2(count)=lambda_2;
            lambda_2=lambda_t-AAlambda;
            resi1=norm(lambda_2,'fro')/norm(lambda_t,'fro');
            lambda_t=lambda_t+z_tplus1-y_tplus1;
    
            resi=norm(abs(y_tplus1)-sqrt(Z),'fro')/norm(sqrt(Z),'fro');
            residual_DR_y(num_im, count)=resi;
            
            
            
            
            residual_DR_lambda(num_im, count)=resi1;
            %fprintf('count=%d\n norm(lambda_2)=%f\n resi=%f\n',count,resi1,resi)
            
            %IM_resi=norm(x_t,'fro')^2+norm_IM^2-2* abs(sum(sum(IM.*conj(x_t))));
            %dev(count)=sqrt(norm(IM_resi))/norm_IM;
            %fprintf('count=%d\n lambda_2=%f\n resi=%f\n', count,lambda_2,resi)
            count=count+1;
            %figure
            %imshow(abs(x_t)/255)
            %pause(0.5)
    
    end
    
 
  num_im=num_im+1  
    
end
figure
[X,Y]=meshgrid(SUB_IM,1:MaxIter);
surf(X,Y,residual_DR_y')
xlabel('size of image n-by-n')
ylabel('iter')
zlabel('residual')
title('gamma=1.5 lambda_0=abs(Y)')


figure

[X,Y]=meshgrid(SUB_IM,1:MaxIter);
surf(X,Y,residual_DR_lambda')
xlabel('size of image n-by-n')
ylabel('iter')
zlabel('||lambda_2||_2/||lambda||_2')
title('gamma=1.5 lambda_0=abs(Y)')




    
    
    