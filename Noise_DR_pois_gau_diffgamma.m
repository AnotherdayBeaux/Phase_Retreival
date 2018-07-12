%%%% Compare  DR_gau, DR_pois with difference gamma




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
MaxIter=500;
Tol=1e-14;

SUBIM=156;

relative_DR_y=zeros(2,MaxIter);
DENSR=0.02:0.04:0.3; % noise level used in the test
Gamma=[5,10,100,1000];
l_g=length(Gamma);
NSR_V=zeros(1,length(DENSR));
end_error=zeros(2*l_g,length(DENSR));


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


lambda_0=Y+rand(size(Y))*4;
%%% add gaussian random noise

trial=1;
for DeNSR=DENSR
    

%tao=(norm(Y,'fro')^2*DeNSR^2)/(size(Y,1)*size(Y,2)); % variance of noise
noise=normrnd(0,1,size(Y,1),size(Y,2));
NSR= norm(noise,'fro')/norm(Y,'fro');
B=abs(Y)+noise*DeNSR/NSR; % new abs(Y)
Z=B.^2; 
%%% verify NSR 
NSR= norm(B-abs(Y),'fro')/norm(Y,'fro');
fprintf('DeNSR=%f,\n NSR=%f\n',DeNSR,NSR);
NSR_V(trial)=NSR;
%%%% end verify   


    gamma_count=1;
    for gamma=Gamma
    %DR_Pois
    lambda_t=lambda_0;
    count_DR_pois=1;   
    resi=1;
    
    while resi > Tol && count_DR_pois < MaxIter
    
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
    
            resi=norm(abs(y_tplus1)-abs(Y),'fro')/norm(abs(Y),'fro');
            
            ee=norm(Y(:)'*AAlambda(:))/(Y(:)'*AAlambda(:));
            rel_y=norm(ee*AAlambda-Y,'fro')/norm(Y,'fro');
            relative_DR_y(1,count_DR_pois)=rel_y;
            %fprintf('count_1=%d\n resi=%f\n rel_y=%f\n', count_DR_pois,resi,rel_y);
            count_DR_pois=count_DR_pois+1;
    
    end
    
    
    
   
        
    %DR_gau
    lambda_t=lambda_0;
    count_DR_gau=1;   
    resi=1;
    while resi > Tol && count_DR_gau < MaxIter
    
            x_t=Nos_ifft(lambda_t,os_rate,num_mask,mask);
            AAlambda=Nos_fft(x_t,os_rate,num_mask,mask);
            y_tplus1=AAlambda;

            ylambda_k=2*y_tplus1-lambda_t;
            Q=abs(ylambda_k);
            % POISSON
            %rot_z=(Q/gamma+sqrt(Q.^2/gamma^2+8*(2+1/gamma)*(abs(Y).^2)))/(4+2/gamma);
            
            % GAUSSIAN
            rot_z=(B+1/gamma*Q)/(1+1/gamma);
            
            ang_z_real=real(ylambda_k)./Q;
            ang_z_imag=imag(ylambda_k)./Q;
            z_tplus1=rot_z.*ang_z_real+1i*rot_z.*ang_z_imag;
    
            lambda_t=lambda_t+z_tplus1-y_tplus1;
    
            resi=norm(abs(y_tplus1)-abs(Y),'fro')/norm(abs(Y),'fro');
            
            ee=norm(Y(:)'*AAlambda(:))/(Y(:)'*AAlambda(:));
            rel_y=norm(ee*AAlambda-Y,'fro')/norm(Y,'fro');
            relative_DR_y(2,count_DR_gau)=rel_y;
            %fprintf('count_3=%d\n resi=%f\n rel_y=%f\n', count_DR_gau,resi,rel_y);
            count_DR_gau=count_DR_gau+1;
    end
           
      end_error(gamma_count,trial)=min(relative_DR_y(2,count_DR_gau-5:count_DR_gau-1));
      end_error(gamma_count+l_g,trial)=min(relative_DR_y(1,count_DR_pois-5:count_DR_pois-1));
      
      gamma_count=gamma_count+1;
        
    end
        
        
        
        
        


trial=trial+1;             
end

figure
plot(NSR_V,end_error(1,:),'rs-')
hold on
plot(NSR_V,end_error(2,:),'gs-')
hold on
plot(NSR_V,end_error(3,:),'bs-')
hold on
plot(NSR_V,end_error(4,:),'cs-')
hold on
%plot(NSR_V,end_error(5,:),'ks-')
%hold on


plot(NSR_V,end_error(5,:),'ro:')
hold on
plot(NSR_V,end_error(6,:),'go:')
hold on
plot(NSR_V,end_error(7,:),'bo:')
hold on
plot(NSR_V,end_error(8,:),'co:')
%hold on
%plot(NSR_V,end_error(10,:),'ko:')
grid on
xlabel('NSR')
ylabel('Rel Error')
legend('DR gau gamma=0.5','DR gau gamma=1','DR gau gamma=1.5','DR gau gamma=2.5','DR gau gamma=4','DR pois gamma=0.5','DR pois gamma=1','DR pois gamma=1.5', 'DR pois gamma=2.5','DR pois gamma=4');   
%Gamma=[0.5,1,1.5,2.5,4];

legend('DR gau gamma=5','DR gau gamma=10','DR gau gamma=100','DR gau gamma=1000','DR pois gamma=5','DR pois gamma=10','DR pois gamma=100', 'DR pois gamma=1000');   
%Gamma=[0.5,1,1.5,2.5,4];


