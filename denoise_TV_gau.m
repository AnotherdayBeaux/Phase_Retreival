%%%% total variation denoising effect. poisson noise case




Bigname='/Users/Beaux/Desktop/phase retrieval'   ; % location of image

str1=[Bigname,'/image_lib/Cameraman.bmp'];   %% consider Im= Cameraman+i* barbera
str2=[Bigname,'/image_lib/Cameraman.bmp'];

os_rate=2;  %% oversampling ratio
TVswitch=0;    % TV   1==ON
               %      0==OFF
               
tv=0.1; % tv * ||grad_image||_1
type_noise = 0;  % =0 pois
                 % =1 gau
Noiseswitch = 1; % 1= add poisson noise to |Y|^2;        
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

MaxIter=3000;
maxiter=4000;
No_Tv=200;
Tv=10;

Tol=1e-14;
toler=1e-8;  % tolerance in total variation denoising
tau= 1/8; % usually < 1/8;

SUBIM=256;

relative_DR_y=zeros(2,MaxIter);
DENSR=20;      %DENSR 1:5
Gamma=5;       %Gamma=[1, 2.5,4];
l_g=length(Gamma);
NSR_V=zeros(1,length(DENSR));
end_error=zeros(2*l_g,length(DENSR));


[IM,Y,Z,SNR,mask]=proc_TV(str1,str2,os_rate,num_mask,add_phase,...
                      sector_alpha,sector_beta,sector_mode,coer,SUBIM);
 if Noiseswitch == 0
     Z=abs(Y).^2;
 end
 
 s
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


%%% add gaussian random noise
trial=1;
for DeNSR=DENSR
    
    
    switch type_noise
        case 1    % gau noise add to b^2
        noise=normrnd(0,abs(Y),size(Y,1),size(Y,2));
        NSR= norm(noise,'fro')/norm(Y,'fro');
        Z=abs(Y).^2+noise*DeNSR^2/NSR; % new abs(Y)
        B=sqrt(abs(Z)); 
        %%% verify NSR 
        NSR= norm(B-abs(Y),'fro')/norm(Y,'fro');
        fprintf('DeNSR=%f,\n NSR=%f\n',DeNSR,NSR);
        NSR_V(trial)=NSR;
        case 0    % pois noise
        [IM,Y,Z,SNR,mask]=proc_TV(str1,str2,os_rate,num_mask,add_phase,...
                      sector_alpha,sector_beta,sector_mode,DeNSR,SUBIM);
        B=sqrt(abs(Z));      
        %%% verify NSR 
        NSR= norm(B-abs(Y),'fro')/norm(Y,'fro');
        fprintf('DeNSR=%f,\n NSR=%f\n',DeNSR,NSR);
            
    end 
        
        
        
    if Noiseswitch == 0
        Z=abs(Y).^2;
    end 

    %%% lambda_0=Y+rand(size(Y))*4;
    lambda_0=rand(size(Y))+1i*rand(size(Y));
    gamma_count=1;
    for gamma=Gamma
        
    lambda_t=lambda_0;
    count_DR=1;   
    rel_y=1;
    
    while rel_y > NSR*0 && count_DR < No_Tv
    
            x_t=Nos_ifft(lambda_t,os_rate,num_mask,mask); 
            AAlambda=Nos_fft(x_t,os_rate,num_mask,mask);
            y_tplus1=AAlambda;

            ylambda_k=2*y_tplus1-lambda_t;
            Q=abs(ylambda_k);
            % POISSON
            %rot_z=(Q/gamma+sqrt(Q.^2/gamma^2+8*(2+1/gamma)*(Z)))/(4+2/gamma);
            
            % GAUSSIAN
            rot_z=(sqrt(abs(Z))+1/gamma*Q)/(1+1/gamma);
            
            ang_z_real=real(ylambda_k)./Q;
            ang_z_imag=imag(ylambda_k)./Q;
            z_tplus1=rot_z.*ang_z_real+1i*rot_z.*ang_z_imag;
    
            lambda_t=lambda_t+z_tplus1-y_tplus1;
   
            
            ee=norm(Y(:)'*AAlambda(:))/(Y(:)'*AAlambda(:));
            rel_y=norm(ee*AAlambda-Y,'fro')/norm(Y,'fro');
            relative_DR_y(1,count_DR)=rel_y;
            fprintf('count_DR=%d\n NSR=%f\n rel_y=%f\n', count_DR,NSR,rel_y);
            count_DR=count_DR+1;

    end
            ee=norm(IM(:)'*x_t(:))/(IM(:)'*x_t(:));
            x_t_=ee*x_t;
            switch type_noise
                case 1   %gau
                figure
                imshow(real(x_t_)/255)
                tit=['gaunotv',',NSR=',num2str(NSR),',rely=',num2str(rel_y)];
                title(tit)
                case 0   %pois
                figure
                imshow(real(x_t_)*DeNSR/255)
                tit=['poisnotv',',NSR=',num2str(NSR),',rely=',num2str(rel_y)];
                title(tit)
            end
            %figure
            %imshow(imag(x_t_)/255)
    %%%%%%%%%
    %%%%%%%%% adding tv denoising after it gets close enough to optimal
        %[x_t,count,err]=update_x_t_TV(x_t,tv,tau,toler, maxiter);       
        %fprintf('count=%d \n err=%s\n', count,err)
        %AAlambda=Nos_fft(x_t,os_rate,num_mask,mask);
        %ee=norm(Y(:)'*AAlambda(:))/(Y(:)'*AAlambda(:));
        %rel_y=norm(ee*AAlambda-Y,'fro')/norm(Y,'fro');
        
        while   count_DR < (No_Tv+Tv)
    
            x_t=Nos_ifft(lambda_t,os_rate,num_mask,mask); 
            %%%
            %%%  solving 1/2 ||y-x_t||^2 + 1/tv * ||y||_BV
            %%%
                
            [x_t,count,err]=update_x_t_TV(x_t,tv,tau,toler, maxiter);
            
            %%%%
            %%%%
            
            
            AAlambda=Nos_fft(x_t,os_rate,num_mask,mask);
            y_tplus1=AAlambda;

            ylambda_k=2*y_tplus1-lambda_t;
            Q=abs(ylambda_k);
            % POISSON
            rot_z=(Q/gamma+sqrt(Q.^2/gamma^2+8*(2+1/gamma)*(Z)))/(4+2/gamma);
            
            % GAUSSIAN
            %rot_z=(sqrt(abs(Z))+1/gamma*Q)/(1+1/gamma);
            
            ang_z_real=real(ylambda_k)./Q;
            ang_z_imag=imag(ylambda_k)./Q;
            z_tplus1=rot_z.*ang_z_real+1i*rot_z.*ang_z_imag;
    
            lambda_t=lambda_t+z_tplus1-y_tplus1;
   
            
            ee=norm(Y(:)'*AAlambda(:))/(Y(:)'*AAlambda(:));
            rel_y=norm(ee*AAlambda-Y,'fro')/norm(Y,'fro');
            relative_DR_y(1,count_DR)=rel_y;
            fprintf('count_DR_tv=%d\n rel_y=%f\n count=%d\n error=%s\n', count_DR,rel_y,count,err);
            count_DR=count_DR+1;

        end
    %%%%%%%%%%%%%
    %%%%%%%%%%%%%

       ee=norm(IM(:)'*x_t(:))/(IM(:)'*x_t(:));
       x_t_tv=ee*x_t;
       switch type_noise
           case 1 % gau
            figure
            imshow(real(x_t_tv)/255)
            tit=['gautv,applied',num2str(Tv),',NSR=',num2str(NSR),',rely=',num2str(rel_y)];
            title(tit)
           case 0 % pois
            figure
            imshow(real(x_t_tv)*DeNSR/255)
            tit=['poistv,applied',num2str(Tv),',NSR=',num2str(NSR),',rely=',num2str(rel_y)];
            title(tit)  
       end
       
       %figure
       %imshow(imag(x_t_tv)/255)
  
    
      
      gamma_count=gamma_count+1;
        
    end
        
        
        
   
trial=trial+1;             
end


noise_IM=x_t_-x_t_tv;
figure
imshow(real(noise_IM))
title('Notv-Tv')