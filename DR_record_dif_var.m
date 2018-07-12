% in this program, we record different variables of DR method. 
% residual error is defined as || |A'Ax|-b||/||b||
% rel error is defined as ||y_tplus1-y_t||/||y_t||
% dev_from_lambda and dev_from_y are of course the same.




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

SUBIM=128;
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
%for i = 1: 10
% DR
%S=abs(Y).*exp(2i*pi*rand(size(Y)));
%z_0=Nos_ifft(S,os_rate,num_mask,mask);
%z_0=Nos_fft(z_0,os_rate,num_mask,mask);
lambda_0=rand(size(Y)).*exp(2i*pi*rand(size(Y)));
%lambda_0=2*z_0-S;
% ADMM
z_0=Nos_ifft(lambda_0,os_rate,num_mask,mask);
z_0=Nos_fft(z_0,os_rate,num_mask,mask);
lambda_00=z_0;
%3
%lambda_0=Y+rand(size(Y))*1e-1;
lambda_2=0;




MaxIter=12000;
Tol=1e-12;
Gamma=2:10:100;    % to be modified
trial=0;

Bigname=[Bigname,'/code/1/DR_test_liguoyi_paper'];
mkdir(Bigname)
for gamma=Gamma
    residual_DR_y=zeros(MaxIter,1);
    residual_DR_lambda=zeros(MaxIter,1);
    rel_DR_y=zeros(MaxIter,1);
    rel_DR_lambda=zeros(MaxIter,1);
    %%%%%%%%%%%%%%
    trial=trial+1;
    %%%%%%%%%%%%%%
    lambda_t=lambda_0;
    count=1;  
    Lambda_2=zeros(MaxIter,1);
    dev_from_lambda=zeros(MaxIter,1);
    dev_from_y=zeros(MaxIter,1);
    y_t=ones(size(lambda_0));
    resi=1;
    resi1=1;
    tic;
 gamma_2nd=10000;
    while resi1 > Tol && count < MaxIter
    
            %%%%%% y_t+1 = argmin 1/2||y-A'A\lambda_t||^2+1/2gamma ||y-lambda_t||^2 %%%%%%%%
            x_t=Nos_ifft(lambda_t,os_rate,num_mask,mask);
            AAlambda=Nos_fft(x_t,os_rate,num_mask,mask);
            %Resi=norm(y_tplus1,'fro')^2+norm(AAlambda,'fro')^2-2*abs(sum(sum(y_tplus1.*conj(AAlambda))));
            %resi1(count)=sqrt(norm(Resi))/norm(y_tplus1,'fro');
            
            
            %y_tplus1=AAlambda; % rho --> 0
            y_tplus1=(gamma*AAlambda+lambda_t)/(1+gamma);
    
            %%%%% z_t+1 = argmin |z|^2-b log(|z|^2)+1/2gamma*||2y_t+1-lambda_t-z||^2
            ylambda_k=2*y_tplus1-lambda_t;
            Q=abs(ylambda_k);
            % POISSON
            rot_z=(Q/gamma_2nd+sqrt(Q.^2/gamma_2nd^2+8*(2+1/gamma_2nd)*Z))/(4+2/gamma_2nd);
            %rot_z=abs(Y);   % when gamma_2 --> inf
            
            % GAUSSIAN
            %rot_z=(sqrt(Z)+1/gamma*Q)/(1+1/gamma);
            ang_z_real=real(ylambda_k)./Q;
            ang_z_imag=imag(ylambda_k)./Q;

            z_tplus1=rot_z.*ang_z_real+1i*rot_z.*ang_z_imag;
    
            %%%% lambda_t+1=lambda_t+z_t+1-y_t+1;
            lambda_2=norm(lambda_t-AAlambda)/norm(lambda_t);
            
            Lambda_2(count)=lambda_2;
            lambda_tpre=lambda_t;
            lambda_t=lambda_t+z_tplus1-y_tplus1;
            
            
            
            %%%%%
            resi=norm(abs(y_tplus1)-sqrt(Z),'fro')/norm(sqrt(Z),'fro');
            residual_DR_y(count)=resi;
            %%%%%
            resi1=norm(abs(lambda_t)-sqrt(Z),'fro')/norm(sqrt(Z),'fro');
            residual_DR_lambda(count)=resi1;
            %%%%%
            
            ee=norm(y_t(:)'*y_tplus1(:))/(y_t(:)'*y_tplus1(:));
            rel_y=norm(ee*y_tplus1-y_t,'fro')/norm(y_t,'fro');
            
            rel_DR_y(count)=rel_y;
            y_t=y_tplus1;
            %%%%%
            
            ee=norm(lambda_tpre(:)'*lambda_t(:))/(lambda_tpre(:)'*lambda_t(:));
            rel_lambda=norm(ee*lambda_t-lambda_tpre,'fro')/norm(lambda_tpre,'fro');
            rel_DR_lambda(count)=rel_lambda;
            %%%%%
            
            ee=norm(IM(:)'*x_t(:))/(IM(:)'*x_t(:));
            IM_resi_lambda=norm(ee*x_t-IM,'fro')/norm_IM;
            dev_from_lambda(count)=IM_resi_lambda;
            %%%%%
            
            y_x_t=Nos_ifft(y_tplus1,os_rate,num_mask,mask);
            ee=norm(IM(:)'*y_x_t(:))/(IM(:)'*y_x_t(:));
            IM_resi_y=norm(ee*y_x_t-IM,'fro')/norm_IM;
            dev_from_y(count)=IM_resi_y;
            
            %fprintf('count=%d\n lambda_2=%f\n resi=%f\n', count,lambda_2,resi)
            count=count+1;
    
    end
    time=toc;
    fixed_point_verify=0;
    if lambda_2 <1e-6
        fixed_point_verify=1;
    end
    
    fprintf('DR \n gamma= %f\n Noise=%d\n number of iter=%d\n Tol=%s\n end_resi=%s \nfixed_point=%d\n time=%f\n'...
                 , gamma, Noiseswitch, (count-1),Tol,resi,fixed_point_verify,time);
    count_DR=count;
    
    figure
    semilogy(1:count_DR-1,residual_DR_y(1:count_DR-1),'g')
    hold on
    semilogy(1:count_DR-1,residual_DR_lambda(1:count_DR-1),'b')
    hold on
    semilogy(1:count_DR-1,Lambda_2(1:count_DR-1),'r')
    hold on 
    semilogy(1:count_DR-1,rel_DR_y(1:count_DR-1),'g-.')
    hold on
    semilogy(1:count_DR-1,rel_DR_lambda(1:count_DR-1),'b-.')
    
    hold on
    semilogy(1:count_DR-1,dev_from_y(1:count_DR-1),'g:')
    hold on
    semilogy(1:count_DR-1,dev_from_lambda(1:count_DR-1),'b:')
    xlabel('iter')
    ylabel('rel/residual error')
    tit=['DR rel/residual error when gamma=',num2str(gamma)];
    legend('res y','res lambda',...
        'rel lambda_2','rel y','rel lambda', 'rel IMAGE in y', 'rel IMAGE in lambda')
    title(tit)
   
    name=[Bigname,'/DR_poisson_gamma=',num2str(gamma*1000)];
    saveas(gcf,name,'png')
    close
    
    
%end

end


