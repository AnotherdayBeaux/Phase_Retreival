
%compare chen's raw DR with chen's+A^*A




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
MaxIter=20000;
first_iter=200;
Tol=1e-14;

SUBIM=56;

relative_DR_y=zeros(2,MaxIter);
DENSR=1:0.5:5; % noise level used in the test
NSR_V=zeros(1,length(DENSR));
end_error=zeros(2,length(DENSR));


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


%lambda_0=Y+rand(size(Y))*0.1;                     % near opt point
lambda_0=rand(size(Y)).*exp(1i*rand(size(Y)))*2;   % random ini guess
%%% add gaussian random noise

trial=1;
for DeNSR=0*DENSR(1:1)
    

%tao=(norm(Y,'fro')^2*DeNSR^2)/(size(Y,1)*size(Y,2)); % variance of noise
noise=normrnd(0,1,size(Y,1),size(Y,2));
NSR= norm(noise,'fro')/norm(Y,'fro');
Z=abs(Y).^2+noise*DeNSR^2/NSR; % new abs(Y)
B=sqrt(Z); 
%%% verify NSR 
NSR= norm(B-abs(Y),'fro')/norm(Y,'fro');
fprintf('DeNSR=%f,\n NSR=%f\n',DeNSR,NSR);
NSR_V(trial)=NSR;
%%%% end verify

    % albert chen DR
    lambda_t=lambda_0;
    count_ini=1;  
    
    resi=1;
    while resi > Tol && count_ini < MaxIter
        
            
            lambda_tt=B.*real(lambda_t)./abs(lambda_t)+1i*B.*imag(lambda_t)./abs(lambda_t);
            x_t=Nos_ifft(2*lambda_tt-lambda_t,os_rate,num_mask,mask);
            AAlambda=Nos_fft(x_t,os_rate,num_mask,mask);
            lambda_t=lambda_t+AAlambda-lambda_tt;
        
            ee=norm(Y(:)'*AAlambda(:))/(Y(:)'*AAlambda(:));
            rel_y=norm(ee*AAlambda-Y,'fro')/norm(Y,'fro');
            
            resi=norm(abs(AAlambda)-abs(Y),'fro')/norm(abs(Y),'fro');
            
            relative_DR_y(1,count_ini)=rel_y;
            %fprintf('count_ini=%d\n resi=%f\n rel_y=%f\n', count_ini,resi,rel_y)
            count_ini=count_ini+1;
    
    end
    
   
    
    
% albert chen with A^*A modify    
lambda_t=lambda_0;    
count_AA=1;
resi=1;
        while resi > Tol && count_AA < MaxIter
            lambda_t=Nos_ifft(lambda_t,os_rate,num_mask,mask);
            lambda_t=Nos_fft(lambda_t,os_rate,num_mask,mask);
            %%%%%% y_t+1 = argmin 1/2||y-A'A\lambda_t||^2+1/2gamma ||y-lambda_t||^2 %%%%%%%%
            lambda_tt=B.*real(lambda_t)./abs(lambda_t)+1i*B.*imag(lambda_t)./abs(lambda_t);
            x_t=Nos_ifft(2*lambda_tt-lambda_t,os_rate,num_mask,mask);
            AAlambda=Nos_fft(x_t,os_rate,num_mask,mask);
            lambda_t=lambda_t+AAlambda-lambda_tt;
        
            ee=norm(Y(:)'*AAlambda(:))/(Y(:)'*AAlambda(:));
            rel_y=norm(ee*AAlambda-Y,'fro')/norm(Y,'fro');
            
            resi=norm(abs(AAlambda)-abs(Y),'fro')/norm(abs(Y),'fro');
            
            relative_DR_y(2,count_AA)=rel_y;
            %fprintf('count_AA=%d\n resi=%f\n rel_y=%f\n', count_AA,resi,rel_y)
            count_AA=count_AA+1;
       
        end
    
figure 
semilogy(1:10000-1, relative_DR_y(1,1:10000-1),'r')
hold on
semilogy(1:10000-1, relative_DR_y(2,1:10000-1),'r--')
tit=['NSR=',num2str(NSR)];
title(tit);
xlabel('iter')
ylabel('rel AAlambda')
legend('chen','mod chen')
Name=['/Users/Beaux/Desktop/phase retrieval/code/2/compare_chen_mchen/',tit];
mkdir(Name)
Name_fig=[Name,'/fig'];
saveas(gcf,Name_fig,'png')
close   
end_errchen=min(relative_DR_y(1,count_ini-5:count_ini-1));
end_error(1,trial)=end_errchen;
end_errmod_chen=min(relative_DR_y(2,count_AA-5:count_AA-1));
end_error(2,trial)=end_errmod_chen;
conf=[Name,'/configration.txt'];
        fid=fopen(conf,'wt');
        fprintf(fid, 'imagesize=%d\n NSR=%f\n end_errorchen=%f\n end_err_modchen=%f\n'...
            ,SUBIM,NSR,end_errchen,end_errmod_chen );
        fclose(fid);

trial=trial+1;             
end

figure
plot(NSR_V,end_error(1,:),'rs-')
hold on
plot(NSR_V,end_error(2,:),'gs-')
grid on
xlabel('NSR')
ylabel('Rel Error')
legend('ChenDR','Mod ChenDR')





