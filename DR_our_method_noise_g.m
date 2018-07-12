







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
Tol=1e-14;
SUB_IM=150;


gamma=100;
rho=1/gamma;

residual_DR_y=zeros(1,MaxIter);
rel_DR_y=zeros(1,MaxIter);
dis_DR_Y = zeros(1,MaxIter);
dis_DR_Y_n=zeros(1,MaxIter);





[IM,Y,Z,SNR,mask]=proc_(str,os_rate,num_mask,add_phase,...
                      sector_alpha,sector_beta,sector_mode,coer,SUB_IM,1);
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

% DR initial
%lambda_0=Y+rand(size(Y))*0.1;
lambda_0=rand(size(Y)).*exp(1i*rand(size(Y)))*2;           % random ini guess
%lambda_0=rand(size(Y)).*Y.*exp(1i*1e-2*rand(size(Y)))*2;  % near opt phase

%%% add gaussian noise to Y %%%

NSR=0.0;
a_r=NSR*norm(Y)/sqrt(size(Y,1)*size(Y,2));

B=abs(Y)+a_r*10*randn(size(Y));
%%% verify NSR 
NSR_V= norm(B-abs(Y))/norm(Y);
fprintf('NSR=%f,\n NSR_V=%f\n',NSR,NSR_V);
pause(1)



    lambda_t=lambda_0;
    count_ini=1;  
    
    resi=1;
    y_t=zeros(size(lambda_0));
    while resi > Tol && count_ini < 20000
           
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
            rot_z=(B+1/gamma*Q)/(1+1/gamma);
            ang_z_real=real(ylambda_k)./Q;
            ang_z_imag=imag(ylambda_k)./Q;

            z_tplus1=rot_z.*ang_z_real+1i*rot_z.*ang_z_imag;
    
            %%%% lambda_t+1=lambda_t+z_t+1-y_t+1;
            %lambda_2=norm(abs(y_tplus1)-abs(lambda_t),'fro')/norm(abs(y_tplus1),'fro');
            
            %Lambda_2(count)=lambda_2;
            %lambda_2=lambda_t-AAlambda;
            %resi1=norm(lambda_2,'fro')/norm(lambda_t,'fro');
            lambda_tplus1=lambda_t+z_tplus1-y_tplus1;
    
            resi=norm(abs(y_tplus1)-abs(Y),'fro')/norm(abs(Y),'fro');
            residual_DR_y(count_ini)=resi;
            %%ee=norm(y_tplus1(:)'*Y(:))/(y_tplus1(:)'*Y(:));
            %%err=norm(ee*Y-y_tplus1,'fro')/norm(Y,'fro');
            dist=sqrt(1-norm(y_tplus1(:)'*Y(:),'fro')^2/norm(y_tplus1,'fro')^2/norm(Y,'fro')^2);
            dis_DR_Y_n(count_ini)=dist;
            
            
            %%dis_DR_Y(count_ini)=err;
            
            
            %ee=norm(y_tplus1(:)'*B(:))/(y_tplus1(:)'*B(:));
            %err2=norm(ee*B-y_tplus1,'fro')/norm(B,'fro');
            %dis_DR_Y_n(count_ini)=err2;
            
            lambda_t=lambda_tplus1;
            
            
            ee=norm(y_t(:)'*y_tplus1(:))/(y_t(:)'*y_tplus1(:));
            rel_y=norm(ee*y_tplus1-y_t,'fro')/norm(y_t,'fro');
            rel_DR_y(count_ini)=rel_y;
            y_t=y_tplus1;
            
            
            
            
            
            %residual_DR_lambda(num_im, count)=resi1;
            %fprintf('count=%d\n norm(lambda_2)=%f\n resi=%f\n',count,resi1,resi)
            
            %IM_resi=norm(x_t,'fro')^2+norm_IM^2-2* abs(sum(sum(IM.*conj(x_t))));
            %dev(count)=sqrt(norm(IM_resi))/norm_IM;
            fprintf('count_ini=%d\n NSR_V=%f\n resi=%f\n err=%f\n rel_y=%f\n lambda_2=%f\n dist=%f\n'...
                , count_ini,NSR_V,resi,err,rel_y,lambda_2,dist)
            count_ini=count_ini+1;
            %figure
            %imshow(abs(x_t)/255)
            %pause(0.5)
    
     
        
    end

    
    
    
    
    
    
    semilogy(1:count_ini-1,residual_DR_y(1:count_ini-1),'r')
    hold on
    semilogy(1:count_ini-1,dis_DR_Y_n(1:count_ini-1),'g')
    hold on
    %semilogy(1:count_ini-1,dis_DR_Y_n(1:count_ini-1))
    %hold on
    semilogy(1:count_ini-1,NSR_V*ones(1,count_ini-1),'g--')
    legend('residual','disDRY')
    %figure
    %imshow(abs(x_t)/255)
    
    
    
    
    
    