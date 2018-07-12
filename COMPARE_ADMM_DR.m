
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
            
print_im=0; 
ADMM=0;   % 0 admm not run; 1 admm run
% number of masks           
num_of_masks=floor(num_mask/2)+1;
                         
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

SUBIM=256;
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
lambda_0=rand(size(Y)).*exp(2i*pi*rand(size(Y)));
% ADMM
z_0=Nos_ifft(lambda_0,os_rate,num_mask,mask);
z_0=Nos_fft(z_0,os_rate,num_mask,mask);
q_0=2*z_0-lambda_0;
abs_q_0=abs(q_0);

lambda_00
%3
%lambda_0=Y+rand(size(Y))*1e-1;
lambda_2=0;




MaxIter=50000;
Tol=1e-10;
Gamma=0.3:0.2:30;                                     %%%% to be modified
trial=0;
compare_iter_DR=zeros(length(Gamma));
compare_iter_ADMM=zeros(length(Gamma));
 Bigname=[Bigname,'/code/1'];
 mkdir([Bigname,'/ADMM'])
 mkdir([Bigname,'/DR'])
for gamma=Gamma
    residual_DR=zeros(MaxIter);
    residual_ADMM=zeros(MaxIter);
    %%%%%%%%%%%%%%
    trial=trial+1;
    %%%%%%%%%%%%%%
    lambda_t=lambda_0;
    count=1;  
    Lambda_2=zeros(MaxIter,1);
    dev=zeros(MaxIter,1);
    
    resi=1;
    tic;
 %gamma=1.2;
    while resi > Tol && count < MaxIter
    
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
            lambda_2=norm(y_tplus1-lambda_t)/norm(y_tplus1);
            
            Lambda_2(count)=lambda_2;
            lambda_t=lambda_t+z_tplus1-y_tplus1;
    
            resi=norm(abs(y_tplus1)-sqrt(Z),'fro')/norm(sqrt(Z),'fro');
            residual_DR(count)=resi;
            
            IM_resi=norm(x_t,'fro')^2+norm_IM^2-2* abs(sum(sum(IM.*conj(x_t))));
            dev(count)=sqrt(norm(IM_resi))/norm_IM;
            %fprintf('count=%d\n lambda_2=%f\n resi=%f\n', count,lambda_2,resi)
            count=count+1;
    
    end
    time=toc;
    fix_point_verify=0;
    if lambda_2 <1e-6
        fixed_point_verify=1;
    end
    fprintf('DR \n gamma= %f\n Noise=%d\n number of iter=%d\n Tol=%s\n end_resi=%s\nfixed_point=%d\n time=%f\n'...
                 , gamma, Noiseswitch, (count-1),Tol,resi,fixed_point_verify,time);
    compare_iter_DR(trial)=count-1;
    Name=[Bigname,...
            '/DR/add_phase=',num2str(add_phase),'_gamma=',...
            num2str(gamma),...
            '_Noiseswitch=',num2str(Noiseswitch)];

        mkdir(Name)
        conf=[Name,'/configration.txt'];
        fid=fopen(conf,'wt');
        fprintf(fid, 'DR \n gamma= %f\n Noise=%d\n number of iter=%d\n Tol=%s\n end_resi=%s\n,fixed_point=%d\n'...
                 , gamma, Noiseswitch, (count-1),Tol,resi,fixed_point_verify);
        fclose(fid);
    semi=1;
    image_print(Name,count,semi,residual_DR,dev,Lambda_2)
    count_DR=count;
    
    % ADMM
   
    z_00=1/(gamma+1)*q_0+(gamma/(1+gamma)*b).*(real(q_0)./abs_q_0+1i*imag(q_0)./abs_q_0);
    lambda_00=1/2*(abs(z_00)-b).*conj(z_00)./abs(z_00);
    lambda_t=lambda_00;
    z_t=z_00;
    
    
    count=1;
    dev=zeros(MaxIter,1);
    resi=1;
    tic;
    while resi > Tol && count < MaxIter
    
            %%%%%% y_t+1 = argmin 1/2||y-A'A\lambda_t||^2+1/2gamma ||y-lambda_t||^2 %%%%%%%%
            x_t=Nos_ifft(z_t-gamma*lambda_t,os_rate,num_mask,mask);
            AAlambda=Nos_fft(x_t,os_rate,num_mask,mask);
            %Resi=norm(y_tplus1,'fro')^2+norm(AAlambda,'fro')^2-2*abs(sum(sum(y_tplus1.*conj(AAlambda))));
            %resi1(count)=sqrt(norm(Resi))/norm(y_tplus1,'fro');
            y_tplus1=AAlambda;
    
            %%%%% z_t+1 = argmin |z|^2-b log(|z|^2)+1/2gamma*||2y_t+1-lambda_t-z||^2
            ylambda_k=y_tplus1+gamma*lambda_t;
            Q=abs(ylambda_k);
            rot_z=(Q/gamma+ sqrt(Q.^2/gamma^2+8*(2+1/gamma)*Z))/(4+2/gamma);
            ang_z_real=real(ylambda_k)./Q;
            ang_z_imag=imag(ylambda_k)./Q;

            z_t=rot_z.*ang_z_real+1i*rot_z.*ang_z_imag;
            %%%% lambda_t+1=lambda_t+z_t+1-y_t+1;
            
            lambda_t=lambda_t+(y_tplus1-z_t)/gamma;
            resi=norm(abs(z_t)-sqrt(Z),'fro')/norm(sqrt(Z),'fro');
            residual_ADMM(count)=resi;
            x_0=Nos_ifft(z_t,os_rate,num_mask,mask);
            IM_resi=norm(x_0,'fro')^2+norm_IM^2-2* abs(sum(sum(IM.*conj(x_0))));
            dev(count)=sqrt(norm(IM_resi))/norm_IM;
            count=count+1;
    end
    time=toc;
    
    fprintf('\n ADMM \n gamma= %f\n Noise=%d\n number of iter=%d\n Tol=%s\n end_resi=%s\n time=%f\n'...
                 , gamma, Noiseswitch, (count-1),Tol,resi,time);
    compare_iter_ADMM(trial)=count-1;
    Name=[Bigname,...
            '/ADMM/add_phase=',num2str(add_phase),'_gamma=',...
            num2str(gamma),...
            '_Noiseswitch=',num2str(Noiseswitch)];

        mkdir(Name)
        conf=[Name,'/configration.txt'];
        fid=fopen(conf,'wt');
        fprintf(fid, 'DR \n gamma= %f\n Noise=%d\n number of iter=%d\n Tol=%s\n end_resi=%s\n'...
                 , gamma, Noiseswitch, (count-1),Tol,resi);
        fclose(fid);
    semi=1;
    image_print(Name,count,semi,residual_ADMM,dev)
     
    
    
    %%%%% semilogy admm/dr
    figure
    semilogy(1:count_DR-1,residual_DR(1:count_DR-1),'g')
    hold on
    semilogy(1:count-1,residual_ADMM(1:count-1),'r')
    title(['gamma= ',num2str(gamma)])
    xlabel('Iter')
    ylabel('|||A^{*}Ax|-b||_2/||b||_2')
    name=[Bigname,'/admm_dr_same_gamma=',num2str(gamma)];
    legend('DR','ADMM')
    saveas(gcf,name,'png')
    close
    %%%%%%%%%%%%%%%%%%%%%%%

end


%%%%% compare different gamma in the DR/ADMM method
figure
plot(Gamma,compare_iter_ADMM,'r')
hold on 
plot(Gamma,compare_iter_DR,'g')
title('iter of admm/DR')
xlabel('gamma')
ylabel('Iter to Tol')
legend('ADMM','DR')
name=[Bigname,'/gamma_verse_iter'];
legend('DR','ADMM')
saveas(gcf,name,'png')
close

