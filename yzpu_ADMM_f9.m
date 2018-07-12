%%%   Id(y)+|z|^2-b log(|z|^2)+alpha||p||_1 + ADMM
%%%   sub to y=z; p=grad u; y= A^*u 
%%%






% pois likelihood 

% gaus likelihood 





Bigname='/Users/Beaux/Desktop/phase retrieval'   ; % location of image


str1=[Bigname,'/image_lib/Cameraman.bmp'];   %% consider Im= Cameraman+i* barbara
str2=[Bigname,'/image_lib/barbara.png'];

os_rate=4;  %% oversampling ratio = 4
rs=1; % rs=downsampling ratio.
% the total sample points collected = rs*os_rate^2

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
%maxiter=2000;
Tol=1e-10;


% tv denoising parameter. 
%toler=1e-12;  % tolerance in total variation denoising
%tau= 1/10;    % usually < 1/8;
%tv=0.1;

SUBIM=156;
type_im=0;  % 0 real
            % 1 str1+i*str1
            % 2 str1+i*str2


relative_DR_x=zeros(1,MaxIter);
DENSR=1;      %DENSR=15:20; % noise level used in the test
l_d=length(DENSR);

Gamma=[1.3, 2.5,4];       %Gamma=[1, 2.5,4];
l_g=length(Gamma);
NSR_V=zeros(1,length(DENSR));

Rec_err=zeros(l_d,l_g);   % reconstruction error

[IM,Y,Z,SNR,mask]=proc_(str1,os_rate,num_mask,add_phase,...
    sector_alpha,sector_beta,sector_mode,coer,SUBIM);
 if Noiseswitch == 0
     Z=abs(Y).^2;
 end
 
 
[Na,Nb]=size(IM);
norm_IM=norm(IM,'fro');
n=Na*Nb; %% dim of data vector
[Y_Na,Y_Nb]=size(Y); %%% This Y_Nb= Nb*os_rate*'num_mask'


RANS=rand(1,Y_Na*Y_Nb);
pos_1=find(RANS<rs);
RANS=zeros(1,Y_Na*Y_Nb);
RANS(pos_1)=1;
RanS=reshape(RANS,Y_Na,Y_Nb);
NRanS=ones(size(RanS));
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
w1w2Matrix=16*((sin(pi*(0:Na-1)/Na)).^2)'*(sin(pi*(0:Nb-1)/Nb)).^2;



% initialize  lambda1_k lambda2_k lambda3_
%             rho1    rho2    rho3
%             y_k     p_k     z_k
%             u_k
%
rho1=1;
rho2=10;
rho3=3000;

lambda1_0=rand(size(Y)).*exp(1i*rand(size(Y))*pi)*0;
lambda2_0=rand([size(IM),2]).*exp(1i*rand([size(IM),2])*pi)*0;
lambda3_0=rand(size(Y)).*exp(1i*rand(size(Y))*pi)*0;

u_0=rand(size(IM)).*exp(1i*rand(size(IM))*pi);
y_0=Nos_fft(u_0,os_rate,num_mask,mask);
z_0=y_0;
p_0=grad_p(u_0);


alpha=200; 
thresh=alpha/rho2; % update p





Bigname='/Users/Beaux/Desktop/phase retrieval/code/8/rand_sample';
mkdir(Bigname);
trial=1;

%%% add poisson random noise
for DeNSR=DENSR
    [IM,Y,Z,SNR,mask]=initial_im(str1,str2,os_rate,num_mask,DeNSR,SUBIM,type_im,RanS);
    if Noiseswitch == 0
        Z=abs(Y).^2;
    end
B=sqrt(abs(Z));  %%% noise case
              
%%% verify NSR 
NSR= norm(B-abs(Y),'fro')/norm(Y,'fro');
fprintf('DeNSR=%f,\n NSR=%f\n',DeNSR,NSR);
NSR_V(trial)=NSR;
%%%% end verify   
 
    %%%%% FOR GAMMA
    name=[Bigname,'/coer=',num2str(DeNSR)];   %% modify 'DeNSR' to 'NSR' 
    % when noiseswitch is changed to 1.
    mkdir(name)
    
 %for gamma=1.2
        %DR_Pois
        lambda1_k=lambda1_0;
        lambda2_k=lambda2_0;
        lambda3_k=lambda3_0;
        
        y_k=y_0;
        z_k=z_0;
        u_k=u_0;
        p_k=p_0;
        
        
        count_DR=1;   
        rel_u=1;
        while rel_u > Tol && count_DR < MaxIter
            
            %%%% 
            % update u
            
             combo_y_lambda3=y_k+lambda3_k/rho3;
             combo_p_lambda2=p_k+lambda2_k/rho2;
             Aylambda3=Nos_ifft(combo_y_lambda3,os_rate,num_mask,mask);
             divplambda2=div_p(combo_p_lambda2);
             
             u_kplus1=ifft2(fft2(rho3/rho2*Aylambda3-divplambda2)./(rho3/rho2+w1w2Matrix));
             
             
            % update y
            inter=rho1*z_k-lambda1_k-lambda3_k;
            Ainter=Nos_ifft(inter,os_rate,num_mask,mask);
            y_kplus1=1/(rho1+rho3)*(Nos_fft(Ainter,os_rate,num_mask,mask)+...
                    rho3*Nos_fft(u_kplus1,os_rate,num_mask,mask));
            
            
            lambda3_kplus1=lambda3_k+rho3*(y_kplus1-Nos_fft(u_kplus1,os_rate,num_mask,mask)); % update lambda3
            
            
            % update p
            grad_u=grad_p(u_kplus1);
            inter=grad_u-lambda2_k/rho2;
            p_kplus1=max(0,abs(inter)-thresh).*inter./abs(inter);
            
            
            
            lambda2_kplus1=lambda2_k+rho2*(p_kplus1-grad_u);   % update lambda2

            % update z
            combo_y_lambda1=y_kplus1+lambda1_k/rho1;
            Q=abs(combo_y_lambda1);
            % POISSON
            rot_z=rho1/(4+2*rho1)*(Q +sqrt(Q.^2+8*Z*(2+rho1)/rho1^2));
            
            % GAUSSIAN
            %rot_z= (sqrt(Z)+rho1*Q)/(1+rho1);
            
            z_kplus1=rot_z.*combo_y_lambda1./Q;
            z_kplus1(isnan(z_kplus1))=0;
            z_kplus1=z_kplus1.*RanS+abs(RanS-1).*combo_y_lambda1;
            
            lambda1_kplus1=lambda1_k+rho1*(y_kplus1-z_kplus1);
   
            
            
            ee=norm(IM(:)'*u_kplus1(:))/(IM(:)'*u_kplus1(:));
            rel_u=norm(ee*u_kplus1-IM,'fro')/norm(IM,'fro');
            relative_DR_x(count_DR)=rel_u;
            fprintf('count_DR=%d\n NSR=%f\n rel_u=%s\n', count_DR,NSR,rel_u);
            
            count_DR=count_DR+1;
            
            lambda1_k=lambda1_kplus1;
            lambda2_k=lambda2_kplus1;
            lambda3_k=lambda3_kplus1;
        
            y_k=y_kplus1;
            z_k=z_kplus1;
            u_k=u_kplus1;
            p_k=p_kplus1;
        end
        ee=norm(IM(:)'*u_k(:))/(IM(:)'*u_k(:));
        u_k_=ee*u_k;
        figure
        imshow(abs(u_k_)*DeNSR/255)
            %imshow(real(x_t_)/255)
        figure
        semilogy(1:count_DR-1, relative_DR_x(1:count_DR-1),'g')
            
        name_=[name,'\',num2str(NSR)];
            %save(gcf,name_,'png')
            %close
        trial=trial+1;
        
        
 end



%[X,Y]=meshgrid(DENSR,Gamma);
%figure
%mesh(X,Y,Rec_err)
%xlabel('Coer')
%ylabel('Gamma')
%zlabel('Error')
%title('TV=0.3')