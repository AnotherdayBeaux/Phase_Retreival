%compare convergence rate of 50-by-50 image and n_of_iteration
% and estimate the lambda_2 of B(G) via variation.
%gaussian 
%poisson






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
Tol=1e-14;
SUB_IM=128;
last=20;  %get convergence rate  2*last interval calculation.


Gamma=[0.2:0.05:1,1.1:0.1:10,10:20];

rel_y_pois=zeros(1,MaxIter);
rel_y_gau =zeros(1,MaxIter);

num_iter_pois=zeros(1,length(Gamma));
num_iter_gau=zeros(1,length(Gamma));

CR_pois=zeros(1,length(Gamma));
CR_gau =zeros(1,length(Gamma));





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
lambda_0=Y+rand(size(Y))*0.1;
%lambda_0=rand(size(Y)).*exp(1i*rand(size(Y)))*2;           % random ini guess
%lambda_0=rand(size(Y)).*Y.*exp(1i*1e-2*rand(size(Y)))*2;  % near opt phase

%%% add gaussian noise to Y %%%

num_gamma=0; % counting the ith gamma 
for gamma=Gamma
    num_gamma=num_gamma+1;
    
    
    
%%%%% gaussian likelihood.
    lambda_t=lambda_0;
    count_ini=1;  
    
    resi=1;
    
    while resi > Tol && count_ini < MaxIter
           
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
            rot_z=(abs(Y)+1/gamma*Q)/(1+1/gamma);
            ang_z_real=real(ylambda_k)./Q;
            ang_z_imag=imag(ylambda_k)./Q;

            z_tplus1=rot_z.*ang_z_real+1i*rot_z.*ang_z_imag;
            lambda_t=lambda_t+z_tplus1-y_tplus1;
    
            resi=norm(abs(y_tplus1)-abs(Y),'fro')/norm(abs(Y),'fro');
            ee=norm(y_tplus1(:)'*Y(:))/(y_tplus1(:)'*Y(:));
            err=norm(ee*Y-y_tplus1,'fro')/norm(Y,'fro');
            rel_y_gau(count_ini)=err;
            %fprintf('count=%d,\n resi=%f,\n',count_ini,resi);
            count_ini=count_ini+1;
          
    end
    figure
    semilogy(1:count_ini-1,rel_y_gau(1:count_ini-1),'r')
    close
    
    cr_g=get_cr(log(rel_y_gau(1:count_ini-1)),last);
    CR_gau(num_gamma)=cr_g;
    num_iter_gau(num_gamma)=count_ini;
    
%%%%%%  poisson likelihood
    lambda_t=lambda_0;
    count_ini=1;  
    
    resi=1;
    
    while resi > Tol && count_ini < MaxIter
           
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
            rot_z=(Q/gamma+ sqrt(Q.^2/gamma^2+8*(2+1/gamma)*Z))/(4+2/gamma);
            %Gaussian likelihood
            %rot_z=(B+1/gamma*Q)/(1+1/gamma);
            ang_z_real=real(ylambda_k)./Q;
            ang_z_imag=imag(ylambda_k)./Q;

            z_tplus1=rot_z.*ang_z_real+1i*rot_z.*ang_z_imag;
            lambda_t=lambda_t+z_tplus1-y_tplus1;
    
            resi=norm(abs(y_tplus1)-abs(Y),'fro')/norm(abs(Y),'fro');
            ee=norm(y_tplus1(:)'*Y(:))/(y_tplus1(:)'*Y(:));
            err=norm(ee*Y-y_tplus1,'fro')/norm(Y,'fro');
            rel_y_pois(count_ini)=err;
            count_ini=count_ini+1;
          
    end
    figure
    semilogy(1:count_ini-1,rel_y_pois(1:count_ini-1),'r')
    close
    
    cr_p=get_cr(log(rel_y_pois(1:count_ini-1)),last);
    CR_pois(num_gamma)=cr_p;
    num_iter_pois(num_gamma)=count_ini;

    
    
    
    fprintf('cr_g=%f,\n cr_p=%f,\n',cr_g,cr_p)
end

b_i= rand(size(Y));
phase_Y=Y./abs(Y);
Y_r=abs(Y);



Nos_ifft(lambda_t,os_rate,num_mask,mask);
            AAlambda=Nos_fft(x_t,os_rate,num_mask,mask);
            
resid=1;           
            
            
while resid > 0.0062
    
    
    resid=Y_r(:)'*b_i(:);
    resid=resid/norm(Y_r,'fro')/norm(b_i,'fro');
    b_i=b_i- resid*Y_r;
    b_i=b_i/norm(b_i,'fro');
    x_0=Nos_ifft(phase_Y.*b_i,os_rate,num_mask,mask);
    c_AAb_i= conj(phase_Y).*Nos_fft(conj(x_0),os_rate,num_mask,mask);
    r_AAb_i=real(c_AAb_i);
    b_i=r_AAb_i;
    
    
    
    
end
lambda_2=norm(b_i,'fro');
            



Rho=1./Gamma;
Rho=flip(Rho);
CR_pois_rho=flip(CR_pois);
CR_gau_rho=flip(CR_gau);
figure
plot(Rho,CR_pois_rho,'g')
hold on
plot(Rho,CR_gau_rho,'r')
xlabel('rho')
ylabel('conv rate')
legend('pois','gau')
tit=['calculated lambda=', num2str(lambda_2)];
title(tit)
Name=['/Users/Beaux/Desktop/phase retrieval/code/3/test the covergence rate/',num2str(SUB_IM),'by',num2str(SUB_IM)];
mkdir(Name)
conf=[Name,'/configration.txt'];
        fid=fopen(conf,'wt');
        fprintf(fid, 'imagesize=%d\n Calclambda2=%f\n ',SUB_IM,lambda_2);
        fclose(fid);
name=[Name,'/cgcr'];
saveas(gcf, name, 'png')
close