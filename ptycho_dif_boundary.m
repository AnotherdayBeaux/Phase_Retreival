%%%%%% ptychography

%%%% compare 0 boundary/ periodic boundary.

%%%% boundary condition = sector constraint. 
%%%% x_t = [A lambda_t]_sector.





Bigname='/Users/Beaux/Desktop/phase retrieval'   ; % location of image

str1=[Bigname,'/image_lib/Cameraman.bmp'];   %% consider Im= Cameraman+i* barbera
str2=[Bigname,'/image_lib/Cameraman.bmp'];

os_rate=2;  %% oversampling ratio
TVswitch=0;    % TV   1==ON
            %      0==OFF
            
 
%%%%%%%%%%%
TV=[0.01,0.05,0.1,0.2, 0.3,0.4, 0.5,0.8,1.5,2]; % tv * ||grad_image||_1, look at 'for tv=TV'
l_tv=length(TV);
%%%%%%%%%%%


Noiseswitch = 1; % 1= add poisson noise to |Y|^2;        
num_mask=1; % 0  0 mask / Id mask
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
MaxIter=300;
maxiter=6000;
Tol=1e-14;
toler=1e-15;  % tolerance in total variation denoising
tau= 1/8; % usually < 1/8;

SUBIM=128;   
overlap=0.5;  %%%% 50% overlaping ratio.
type_im=1;
bd_con=0;





relative_DR_x=zeros(1,MaxIter);
DENSR=10:5:40;      %DENSR=15:20; % noise level used in the test
L=flip(2:6); % subim= 2^l measure the size of each small block.
len_NSR=length(DENSR);

%Gamma=[0.7,0.9,1.3,2,4];       %Gamma=[1, 2.5,4];
gamma=0.7;



len_l=length(L);
NSR_V=zeros(len_l,len_NSR);
end_error=zeros(len_l,len_NSR);


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

Bigname='/Users/Beaux/Desktop/phase retrieval/code/16/ptycho_diff_boundary_00';
mkdir(Bigname);

trial_l=1;
for l=L
    
    subim=2^l;
    q=2*SUBIM/subim;
    [IM,Y,Z,NSR,mask]=ini_ptycho(str1,str2,os_rate,num_mask,coer,SUBIM, type_im,subim,overlap,bd_con);
 
 
    [Na,Nb]=size(IM);
    [Y_Na,Y_Nb]=size(Y);
    
    name_l=[Bigname,'/q=',num2str(q)];
    mkdir(name_l)
    
    

trial_NSR=1;
for DeNSR=DENSR
    [IM,Y,Z,NSR,mask,nor_ptycho]=ini_ptycho(str1,str2,os_rate,num_mask...
        ,DeNSR,SUBIM, type_im,subim,overlap,bd_con);
    R_IM=IM(subim/2+1:Na-subim/2,subim/2+1:Nb-subim/2);   % trucate the genuine image
if Noiseswitch == 0
     Z=abs(Y).^2;
end
B=sqrt(abs(Z));  %%% noise case
              
%%% get& record SNR_intensity DB
%NSR= (norm(B-abs(Y),'fro')/norm(Y,'fro'));
%NSR=-20 *log(SNR)/log(10);
fprintf('DeNSR=%f,\n NSR=%f\n',DeNSR,NSR);
NSR_V(trial_l,trial_NSR)=NSR;
%%%% end verify   

    %%% lambda_0=Y+rand(size(Y))*4;
    lambda_0=rand(size(Y))+1i*rand(size(Y));
    
    
    
    name_NSR=[name_l,'/NSR',num2str(NSR)];
    mkdir(name_NSR)
    %DR_Pois
    lambda_t=lambda_0;
    count_DR=1;   
    rel_x=1;
    
    
    
    switch bd_con
        
        
        % 0 boundary 
        case 0 
            
        while rel_x > Tol && count_DR< MaxIter
    
            x_t=ptycho_ifft(lambda_t,os_rate,num_mask,mask,subim,overlap); 
            x_t=x_t./(nor_ptycho).^2;
            
            
            x__t=zeros(Na,Nb);
            x__t(subim/2+1:Na-subim/2,subim/2+1:Nb-subim/2)=x_t(subim/2+1:Na-subim/2,subim/2+1:Nb-subim/2);
           
            AAlambda=ptycho_fft(x__t,os_rate,num_mask,mask,subim,overlap);
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
   
            R_x_t=x__t(subim/2+1:Na-subim/2,subim/2+1:Nb-subim/2);
            ee=norm(R_IM(:)'*R_x_t(:))/(R_IM(:)'*R_x_t(:));
            rel_x=norm(ee*R_x_t-R_IM,'fro')/norm(R_IM,'fro');
            relative_DR_x(1,count_DR)=rel_x;
            fprintf('count_DR=%d\n NSR=%f\n rel_x=%f\n', count_DR,NSR,rel_x);
            count_DR=count_DR+1;

        end
    
    %%%% periodic boundary 
        case 1
          while rel_x > Tol && count_DR< MaxIter
    
            x_t=ptycho_ifft(lambda_t,os_rate,num_mask,mask,subim,overlap); 
            x__t=x_t./(nor_ptycho.^2);
        
            
            
            x__t=proj_periodic(x__t,subim);
            AAlambda=ptycho_fft(x__t,os_rate,num_mask,mask,subim,overlap);
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
   
            R_x_t=x__t(subim/2+1:Na-subim/2,subim/2+1:Nb-subim/2);
            ee=norm(R_IM(:)'*R_x_t(:))/(R_IM(:)'*R_x_t(:));
            rel_x=norm(ee*R_x_t-R_IM,'fro')/norm(R_IM,'fro');
            relative_DR_x(1,count_DR)=rel_x;
            fprintf('count_DR=%d\n NSR=%f\n rel_x=%f\n', count_DR,NSR,rel_x);
            count_DR=count_DR+1;

         end
    
    
    
   
    end
    
    figure
    semilogy(1:count_DR-1,relative_DR_x(1:count_DR-1),'g')
    xlabel('iter')
    ylabel('relError')
    Title=['NSR=',num2str(NSR),',q=',num2str(q)];
    title(Title)
    
    name_plot=[name_NSR,'/relError'];
    saveas(gcf,name_plot,'png')
    close
 
    ee=norm(R_IM(:)'*R_x_t(:))/(R_IM(:)'*R_x_t(:));
    x_t_=ee*R_x_t;
    figure
    imshow(real(x_t_)*DeNSR/255)
    name_plot=[name_NSR,'/recon_im'];
    saveas(gcf,name_plot,'png')
    close
    end_error(trial_l,trial_NSR)=min(relative_DR_x(count_DR-5:count_DR-1));
    
    
    trial_NSR=trial_NSR+1;   
end

trial_l=trial_l+1;

end
figure
for i = 1: len_l
    
plot(NSR_V(i,:),end_error(i,:),'+-')
hold on

end
grid on
xlabel('NSR')
ylabel('Rec Err')
title('NSR verse Err 0 ')
legend('q=4','q=8','q=16','q=32','q=64')

name_plot=[Bigname,'/NSR verse Err 0'];
saveas(gcf,name_plot,'png')