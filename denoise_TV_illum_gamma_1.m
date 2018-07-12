%%%% total variation denoising effect. 'POISSON' noise case
%%%% in order to verify that when illuminicity is high, the reconstruction
%%%% is good.




Bigname='/Users/Beaux/Desktop/phase retrieval'   ; % location of image

str1=[Bigname,'/image_lib/Cameraman.bmp'];   %% consider Im= Cameraman+i* barbera
str2=[Bigname,'/image_lib/Cameraman.bmp'];

os_rate=2;  %% oversampling ratio

%%%%%%%%%%%
TV=[0.2,0.4,0.6,0.9,1.2,1.5]; % tv * ||grad_image||_1, look at 'for tv=TV'
l_tv=length(TV);
%%%%%%%%%%%


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
            
coer=3;    % regulating the SNR. coer higher, noise higher.
MaxIter=3000;
maxiter=2000;
Tol=1e-14;
toler=1e-12;  % tolerance in total variation denoising
tau= 1/10; % usually < 1/8;

SUBIM=120;

relative_DR_y=zeros(2,MaxIter);
DENSR=[0.3,0.5,0.8,1,1.3];      %DENSR=15:20; % noise level used in the test
l_d=length(DENSR);

RHO=[0.05,1,1.5,3,5,10,15,20];
Gamma=1./RHO;       %Gamma=[1, 2.5,4];
l_g=length(Gamma);
NSR_V=zeros(1,length(DENSR));
end_error=zeros(2*l_g,length(DENSR));

Rec_err=zeros(l_d,l_tv);   % reconstruction error rec_err(nsr,gamma,tv)=err

[IM,Y,Z,SNR,mask]=proc_TV(str1,str2,os_rate,num_mask,add_phase,...
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


%%% add poisson random noise

Bigname='/Users/Beaux/Desktop/phase retrieval/code/7/denoise_tv_illum';
mkdir(Bigname);
trial=1;
for DeNSR=DENSR
    [IM,Y,Z,SNR,mask]=proc_TV(str1,str2,os_rate,num_mask,add_phase,...
                      sector_alpha,sector_beta,sector_mode,DeNSR,SUBIM);
if Noiseswitch == 0
     Z=abs(Y).^2;
end
B=sqrt(abs(Z));  %%% noise case
              
%%% verify NSR 
NSR= norm(B-abs(Y),'fro')/norm(Y,'fro');
fprintf('DeNSR=%f,\n NSR=%f\n',DeNSR,NSR);
NSR_V(trial)=NSR;
%%%% end verify   

    %%% lambda_0=Y+rand(size(Y))*4;
    lambda_0=rand(size(Y))+1i*rand(size(Y));
    gamma_count=1;
    
    
    
    %%%%% FOR GAMMA
    name=[Bigname,'/coer=',num2str(DeNSR)];   %% modify 'DeNSR' to 'NSR' 
    % when noiseswitch is changed to 1.
    mkdir(name)
    
    
    for gamma=Gamma
    %DR_Pois
    lambda_t=lambda_0;
    count_DR=1;   
    rel_y=1;
    while rel_y > 0.1 && count_DR < 20
    
            x_t=Nos_ifft(lambda_t,os_rate,num_mask,mask); 
            AAlambda=Nos_fft(x_t,os_rate,num_mask,mask);
            y_tplus1=AAlambda;

            ylambda_k=2*y_tplus1-lambda_t;
            Q=abs(ylambda_k);
            % POISSON
            rot_z=(Q/gamma+sqrt(Q.^2/gamma^2+8*(2+1/gamma)*(Z)))/(4+2/gamma);
            
            % GAUSSIAN
            %rot_z=(sqrt(Z)+1/gamma*Q)/(1+1/gamma);
            
            ang_z_real=real(ylambda_k)./Q;
            ang_z_imag=imag(ylambda_k)./Q;
            z_tplus1=rot_z.*ang_z_real+1i*rot_z.*ang_z_imag;
    
            lambda_t=lambda_t+z_tplus1-y_tplus1;
   
            
            ee=norm(Y(:)'*AAlambda(:))/(Y(:)'*AAlambda(:));
            rel_y=norm(ee*AAlambda-Y,'fro')/norm(Y,'fro');
            relative_DR_y(1,count_DR)=rel_y;
            fprintf('count_DR=%d\n Coer=%f\n rel_y=%s\n', count_DR,DeNSR,rel_y);
            count_DR=count_DR+1;

    end
            ee=norm(IM(:)'*x_t(:))/(IM(:)'*x_t(:));
            x_t_=ee*x_t;
            figure
            %imshow(real(x_t_)*DeNSR/255)
            imshow(real(x_t_)/255)
            
            name_=[name,'\',num2str(NSR)];
            %save(gcf,name_,'png')
            %close
    %%%%%%%%%
    %%%%%%%%% adding tv denoising after it gets close enough to optimal
     tv_name=[name,'/gam=',num2str(gamma),'_TV',];
     mkdir(tv_name)
        
       tv_count=1;
        for tv=TV
        [x_t,count,err]=update_x_t_TV(x_t,tv,tau,toler, maxiter*iter);   
        AAlambda=Nos_fft(x_t,os_rate,num_mask,mask);
        rel_y=norm(ee*AAlambda-Y,'fro')/norm(Y,'fro');
        fprintf('count=%d \n err=%s\n rel_y=%s\n', count,err,rel_y)
        
        %end
        
        
        
        rel_y=1;
        rel_bef_y=2;
        
        while   count_DR < 80 && rel_bef_y>rel_y
    
            x_t=Nos_ifft(lambda_t,os_rate,num_mask,mask); 
            %%%
            %%%  solving 1/2 ||y-x_t||^2 + 1/tv * ||y||_BV
            %%%
            %%%%
            AAlambda=Nos_fft(x_t,os_rate,num_mask,mask); % get rel_y before 
            %applying update_x_t
            ee=norm(Y(:)'*AAlambda(:))/(Y(:)'*AAlambda(:));
            rel_bef_y=norm(ee*AAlambda-Y,'fro')/norm(Y,'fro');
            fprintf('count_DR_bef_tv=%d\n rel_y=%s\n', count_DR,rel_bef_y);   
            
            
            x_t=update_x_t_TV(x_t,tv,tau,toler, maxiter);

            %%%%
            %%%%
            
            
            AAlambda=Nos_fft(x_t,os_rate,num_mask,mask);
            ee=norm(Y(:)'*AAlambda(:))/(Y(:)'*AAlambda(:));
            rel_y=norm(ee*AAlambda-Y,'fro')/norm(Y,'fro');
            relative_DR_y(1,count_DR)=rel_y;
            fprintf('count_DR_tv=%d\n rel_y=%s\n', count_DR,rel_y);
            y_tplus1=AAlambda;
            
            ylambda_k=2*y_tplus1-lambda_t;
            Q=abs(ylambda_k);
            % POISSON
            %rot_z=(Q/gamma+sqrt(Q.^2/gamma^2+8*(2+1/gamma)*(Z)))/(4+2/gamma);
            
            % GAUSSIAN
            rot_z=(sqrt(Z)+1/gamma*Q)/(1+1/gamma);
            
            ang_z_real=real(ylambda_k)./Q;
            ang_z_imag=imag(ylambda_k)./Q;
            z_tplus1=rot_z.*ang_z_real+1i*rot_z.*ang_z_imag;
    
            lambda_t=lambda_t+z_tplus1-y_tplus1;
   
            
            count_DR=count_DR+1;

        end
    %%%%%%%%%%%%%
    %%%%%%%%%%%%%

       ee=norm(IM(:)'*x_t(:))/(IM(:)'*x_t(:));
       x_t_tv=ee*x_t;
       figure
       imshow(real(x_t_tv)/255)
       name_=[tv_name,'\Re_part_recons'];
       saveas(gcf,name_,'png')
       figure
       imshow(imag(x_t_tv)/255)
       
       name_=[tv_name, '\Im_part_recons'];
       saveas(gcf,name_,'png')
       noise_IM=x_t_-x_t_tv;
       figure
       imshow(real(noise_IM))
       name_=[tv_name,'\error'];
       saveas(gcf,name_,'png')
       close
       close
       close
       
       Rec_err(trial,tv_count)= rel_y;
       tv_count=tv_count+1;
    end
       
     
     
     
     
     
     gamma_count=gamma_count+1;    
    end
        
        
   
trial=trial+1;             
end


[X,Y]=meshgrid(DENSR,TV);
figure
mesh(X,Y,Rec_err')
xlabel('Co')
ylabel('TV')
zlabel('Error')
title('gamma=0.1')