%%% test random sample effect in noiseless case.
%%% 


% pois likelihood 

% gaus likelihood 





Bigname='/Users/Beaux/Desktop/phase retrieval'   ; % location of image


str1=[Bigname,'/image_lib/Cameraman.bmp'];   %% consider Im= Cameraman+i* barbara
str2=[Bigname,'/image_lib/barbara.png'];

os_rate=2;  %% oversampling ratio = 4
rs=0.4; % rs=downsampling ratio.
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

    %%% lambda_0=Y+rand(size(Y))*4;
    lambda_0= 10*rand(size(Y)).*exp(1i*rand(size(Y))*pi);
    gamma_count=1;
    
    
    
    %%%%% FOR GAMMA
    name=[Bigname,'/coer=',num2str(DeNSR)];   %% modify 'DeNSR' to 'NSR' 
    % when noiseswitch is changed to 1.
    mkdir(name)
    
 for gamma=1.2
        %DR_Pois
        lambda_t=lambda_0;
        count_DR=1;   
        rel_x=1;
        while rel_x > Tol && count_DR < MaxIter
    
            x_t=random_sample_ifft(lambda_t,os_rate,num_mask,mask,NRanS); 
            AAlambda=random_sample_fft(x_t,os_rate,num_mask,mask,NRanS);
            y_tplus1=AAlambda;

            ylambda_k=2*y_tplus1-lambda_t;
            Q=abs(ylambda_k);
            % POISSON
            %rot_z=(Q/gamma+sqrt(Q.^2/gamma^2+8*(2+1/gamma)*(Z)))/(4+2/gamma);
            
            % GAUSSIAN
            rot_z=(sqrt(Z)+1/gamma*Q)/(1+1/gamma);
            
            ang_z_real=real(ylambda_k)./Q;
            ang_z_imag=imag(ylambda_k)./Q;
            ang_z_real(isnan(ang_z_real))=0;
            ang_z_imag(isnan(ang_z_imag))=0;
            z_tplus1=rot_z.*ang_z_real+1i*rot_z.*ang_z_imag;
            z_tplus1=z_tplus1.*RanS+abs(RanS-1).*ylambda_k;
            lambda_t=lambda_t+z_tplus1-y_tplus1;
   
            
            
            ee=norm(IM(:)'*x_t(:))/(IM(:)'*x_t(:));
            rel_x=norm(ee*x_t-IM,'fro')/norm(IM,'fro');
            relative_DR_x(count_DR)=rel_x;
            fprintf('count_DR=%d\n NSR=%f\n rel_x=%s\n', count_DR,NSR,rel_x);
            
            
            count_DR=count_DR+1;

        end
            ee=norm(IM(:)'*x_t(:))/(IM(:)'*x_t(:));
            x_t_=ee*x_t;
            figure
            imshow(abs(x_t_)*DeNSR/255)
            %imshow(real(x_t_)/255)
            figure
            semilogy(1:count_DR-1, relative_DR_x(1:count_DR-1),'g')
            
            name_=[name,'\',num2str(NSR)];
            %save(gcf,name_,'png')
            %close
        trial=trial+1;             
 end

end


%[X,Y]=meshgrid(DENSR,Gamma);
%figure
%mesh(X,Y,Rec_err)
%xlabel('Coer')
%ylabel('Gamma')
%zlabel('Error')
%title('TV=0.3')