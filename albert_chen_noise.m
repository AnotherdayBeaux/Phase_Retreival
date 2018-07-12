
%%%% Compare the noise sensitivities between Chen's, PoisDR, GauDR %%%%
%%%%
%%%%

%%%% we explore the two properties: 
%%%% in gaussian noise case.
%%%% 1a. compare chen's and chen's +A^*A. NSR verse accuracy
%%%% 1b. the NSR verse accuracy, with different rhos PoisDR GauDR

%%%% 2ab do the same thing for poisson noise

%%%% 3. os rate effect.
%%%% NSR verse accuracy for different os rate


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
MaxIter=5000;
first_iter=200;
Tol=1e-14;

SUBIM=156;

residual_DR_y=zeros(1,MaxIter);
relative_DR_y=zeros(1,MaxIter);



num_im=1;

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

lambda_0=Y+rand(size(Y))*0.1;
%lambda_0=rand(size(Y)).*exp(1i*rand(size(Y)))*2;           % random ini guess
%lambda_0=rand(size(Y)).*Y.*exp(1i*1e-2*rand(size(Y)))*2;  % near opt phase

%%%
%%% add gaussian random noise
NSR=0.01;
a_r=NSR*norm(Y)/sqrt(size(Y,1)*size(Y,2));

%%% pois noise,



B=abs(Y)+a_r*10*randn(size(Y));
%%% verify NSR 
NSR_V= norm(B-abs(Y))/norm(Y);
fprintf('NSR=%f,\n NSR_V=%f\n',NSR,NSR_V);
pause(2)
%%%% end verify

    lambda_t=lambda_0;
    count_ini=1;  
    
    resi=1;
    AAlambda_pre=zeros(size(lambda_t));
    while resi > Tol && count_ini < MaxIter

            %%%%%% y_t+1 = argmin 1/2||y-A'A\lambda_t||^2+1/2gamma ||y-lambda_t||^2 %%%%%%%%
            lambda_tt=B.*real(lambda_t)./abs(lambda_t)+1i*B.*imag(lambda_t)./abs(lambda_t);
            x_t=Nos_ifft(2*lambda_tt-lambda_t,os_rate,num_mask,mask);
            AAlambda=Nos_fft(x_t,os_rate,num_mask,mask);
            lambda_t=lambda_t+AAlambda-lambda_tt;
        
            ee=norm(AAlambda_pre(:)'*AAlambda(:))/(AAlambda_pre(:)'*AAlambda(:));
            rel_y=norm(ee*AAlambda-AAlambda_pre,'fro')/norm(AAlambda_pre,'fro');
            AAlambda_pre=AAlambda;
            resi=norm(abs(AAlambda)-abs(Y),'fro')/norm(abs(Y),'fro');
            
            
            residual_DR_y(count_ini)=resi;
            relative_DR_y(count_ini)=rel_y;
            fprintf('count_ini=%d\n resi=%f\n rel_y=%f\n', count_ini,resi,rel_y)
            count_ini=count_ini+1;
    
    end
    
    
    
    

    
        while resi > Tol && count_ini < MaxIter
            lambda_t=Nos_ifft(lambda_t,os_rate,num_mask,mask);
            lambda_t=Nos_fft(lambda_t,os_rate,num_mask,mask);
            %%%%%% y_t+1 = argmin 1/2||y-A'A\lambda_t||^2+1/2gamma ||y-lambda_t||^2 %%%%%%%%
            lambda_tt=B.*real(lambda_t)./abs(lambda_t)+1i*B.*imag(lambda_t)./abs(lambda_t);
            x_t=Nos_ifft(2*lambda_tt-lambda_t,os_rate,num_mask,mask);
            AAlambda=Nos_fft(x_t,os_rate,num_mask,mask);
            lambda_t=lambda_t+AAlambda-lambda_tt;
        
            ee=norm(AAlambda_pre(:)'*AAlambda(:))/(AAlambda_pre(:)'*AAlambda(:));
            rel_y=norm(ee*AAlambda-AAlambda_pre,'fro')/norm(AAlambda_pre,'fro');
            AAlambda_pre=AAlambda;
            resi=norm(abs(AAlambda)-abs(Y),'fro')/norm(abs(Y),'fro');
            
            residual_DR_y(count_ini)=resi;
            relative_DR_y(count_ini)=rel_y;
            fprintf('count_ini=%d\n resi=%f\n rel_y=%f\n', count_ini,resi,rel_y)
            count_ini=count_ini+1;
    
       end
    
    
    
    semilogy(1:count_ini-1,residual_DR_y(1:count_ini-1))
    
    
    
    
    
    
   






