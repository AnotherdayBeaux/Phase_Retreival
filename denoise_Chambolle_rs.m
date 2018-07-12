%%%% total variation denoising effect. 'POISSON' noise case

%%%% x_k = A lambda_k;
%%%% x_k = argmin alpha ||x||_1+1/2 ||x-x_k||^2
%%%% y_k = A^*x_k
%%%% z_k = prox(2lambda_k-y_k);
%%%% lambda_k+1= lambda_k+z_k-y_k

%%%% modified: add random sampling feature. 

%%%% We already know the optimal tv(alpha above) lies near 0.3. So we
%%%% assume alpha = 0.3
%%%% in accordance with Chang's ADMM, 3 masks.
%%%%


Bigname='/Users/Beaux/Desktop/phase retrieval'   ; % location of image

str1=[Bigname,'/image_lib/Cameraman.bmp'];   %% consider Im= Cameraman+i* barbera
str2=[Bigname,'/image_lib/Cameraman.bmp'];

os_rate=2;  %% oversampling ratio
TVswitch=0;    % TV   1==ON
            %      0==OFF
%%%%%%%%%%%
tv=0.3; % tv * ||grad_image||_1, look at 'for tv=TV'
%%%%%%%%%%%


Noiseswitch = 1; % 1= add poisson noise to |Y|^2;        
num_mask=4; % 0  0 mask / Id mask
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
            
coer=0.2;    % regulating the SNR. coer higher, noise higher.
MaxIter=10000;
maxiter=8000;
Tol=1e-14;
toler=1e-15;  % tolerance in total variation denoising
tau= 1/8; % usually < 1/8;

SUBIM=50;
type_im=0;

relative_DR_x=zeros(2,MaxIter);
DENSR=20:10:60;      %DENSR=15:20; % noise level used in the test
l_d=length(DENSR);

RHO=[0.5,0.8,1,1.3,10,50,100,600,1000];
len_RHO=length(RHO);
NSR_V=zeros(1,length(DENSR));

% oversampling 
RS=0.7:0.1:1; % random sampling
len_RS=length(RS);
Rec_err=zeros(l_d,len_RHO,len_RS);   % reconstruction error rec_err(SNR_DB,RHO,rs)=err

[IM,Y,Z,SNR,mask]=proc_(str1,os_rate,num_mask,add_phase,...
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

Bigname='/Users/Beaux/Desktop/phase retrieval/code/12/random_sample_chambolle';
mkdir(Bigname);
trial_SNR=1;




for DeNSR=DENSR
  

    trial_RS=1;
    
    for rs=RS
        
        
        RANS=rand(1,Y_Na*Y_Nb);
        pos_1=find(RANS<rs);
        RANS=zeros(1,Y_Na*Y_Nb);
        RANS(pos_1)=1;
        RanS=reshape(RANS,Y_Na,Y_Nb);
        NRanS=ones(size(RanS));
        
        
        [IM,Y,Z,SNR_DB,mask]=initial_im(str1,str2,os_rate,num_mask,DeNSR,SUBIM,type_im,RanS);
        if Noiseswitch == 0
            Z=abs(Y).^2;
        end
        B=sqrt(abs(Z));  %%% noise case
              
        %%% verify NSR 
        NSR= (norm(B-abs(Y),'fro')/norm(Y,'fro'));
        SNR_DB=-10 *log(NSR^2)/log(10);
        %fprintf('DeNSR=%f,\n SNR_DB=%f\n',DeNSR,SNR_DB);



        NSR_V(trial_RS,trial_SNR)=SNR_DB;
        
        name=[Bigname,'/DeNSR=',num2str(DeNSR),',rs=',num2str(rs)];
        mkdir(name)
    
       u_k=Nos_ifft(RanS.*B,os_rate,num_mask,mask);
       lambda_0=Nos_fft(u_k,os_rate,num_mask,mask);
     trial_RHO=1;
     for rho=RHO
        %DR_Pois
        gamma=1/rho;
        lambda_t=lambda_0;
        count_DR=1;   
        rel_x=1;
        
        MaxIter=floor(sqrt(rho)*200);
        while rel_x > NSR*0 && count_DR < MaxIter
    
            x_t=Nos_ifft(RanS.*lambda_t,os_rate,num_mask,mask); 
            AAlambda=Nos_fft(x_t,os_rate,num_mask,mask);
            y_tplus1=AAlambda;

            ylambda_k=2*y_tplus1-lambda_t;
            Q=abs(ylambda_k);
            % POISSON
            rot_z=(Q/gamma+sqrt(Q.^2/gamma^2+8*(2+1/gamma)*(Z)))/(4+2/gamma);
            
            % GAUSSIAN
            %rot_z=(sqrt(Z)+1/gamma*Q)/(1+1/gamma);
            
            z_tplus1=rot_z.*ylambda_k./Q;
            z_tplus1(isnan(z_tplus1))=0;
            z_tplus1=z_tplus1.*RanS;
            z_tplus1=z_tplus1.*RanS+abs(RanS-1).*ylambda_k;
            
            
            lambda_t=lambda_t+z_tplus1-y_tplus1;
          
            
            %ee=norm(Y(:)'*AAlambda(:))/(Y(:)'*AAlambda(:));
            %rel_y=norm(ee*AAlambda-Y,'fro')/norm(Y,'fro');
            %relative_DR_y(1,count_DR)=rel_y;
            %fprintf('count_DR=%d\n NSR=%f\n rel_y=%f\n', count_DR,NSR,rel_y);
            %count_DR=count_DR+1;
            
            
            ee=norm(x_t(:)'*IM(:))/(x_t(:)'*IM(:));
            rel_x=norm(ee*IM-x_t,'fro')/norm(IM,'fro');
            relative_DR_x(1,count_DR)=rel_x;
            rel_xDB=-20*log(rel_x)/log(10);
            %fprintf('count_DR=%d\n SNR_DB=%f\n rel_x=%f\n rel_xDB=%f\n', count_DR,SNR_DB,rel_x,rel_xDB);
            count_DR=count_DR+1;

        end
            ee=norm(IM(:)'*x_t(:))/(IM(:)'*x_t(:));
            x_t_=ee*x_t;
            
            tv_name=[name,'/rho=',num2str(rho)];
            mkdir(tv_name)
            
            
            name_=[name,'\before_TV,SNR_DB=',num2str(SNR_DB),',real part'];
            figure
            imshow(real(x_t_)*DeNSR/255)
            saveas(gcf,name_,'png')
            close
            name_=[name,'\before_TV,SNR_DB=',num2str(SNR_DB),',imag part'];
            figure
            imshow(imag(x_t_)*DeNSR/255)
            saveas(gcf,name_,'png')
            close
            
    %%%%%%%%%
    %%%%%%%%% adding tv denoising after it gets close enough to optimal

        
            [x_t,count,err]=update_x_t_TV(x_t,tv,tau,toler, maxiter);   
            ee=norm(x_t(:)'*IM(:))/(x_t(:)'*IM(:));
            rel_x=norm(ee*IM-x_t,'fro')/norm(IM,'fro');
            %fprintf('count=%d \n err=%s\n', count,err)
        
        while   count_DR < MaxIter+15
    
            x_t=Nos_ifft(RanS.*lambda_t,os_rate,num_mask,mask); 
            %%%
            %%%  solving 1/2 ||y-x_t||^2 + 1/tv * ||y||_BV
            %%%
                
            x_t=update_x_t_TV(x_t,tv,tau,toler, maxiter);

            %%%%
            %%%%
            
            
            AAlambda=Nos_fft(x_t,os_rate,num_mask,mask);
            y_tplus1=AAlambda;

            ylambda_k=2*y_tplus1-lambda_t;
            Q=abs(ylambda_k);
            % POISSON
            %rot_z=(Q/gamma+sqrt(Q.^2/gamma^2+8*(2+1/gamma)*(Z)))/(4+2/gamma);
            
            % GAUSSIAN
            rot_z=(sqrt(Z)+1/gamma*Q)/(1+1/gamma);
            
            z_tplus1=rot_z.*ylambda_k./Q;
            z_tplus1(isnan(z_tplus1))=0;
            z_tplus1=z_tplus1.*RanS;
            z_tplus1=z_tplus1.*RanS+abs(RanS-1).*ylambda_k;
    
            lambda_t=lambda_t+z_tplus1-y_tplus1;
   
            
            ee=norm(IM(:)'*x_t(:))/(IM(:)'*x_t(:));
            rel_x=norm(ee*x_t-IM,'fro')/norm(IM,'fro');
            relative_DR_x(1,count_DR)=rel_x;
            rel_opt=min(relative_DR_x(1,(count_DR-15):count_DR));
            rel_optDB=-20*log(rel_opt)/log(10);
            %fprintf('count_DR_tv=%d\n rel_x=%f\n rel_optDB=%f\n', count_DR,rel_x,rel_optDB);
            count_DR=count_DR+1;

        end
    %%%%%%%%%%%%%
    %%%%%%%%%%%%%
       fprintf('rho=%f\n rs=%f\n SNR_DB=%f\n rel_optDB=%f\n\n',rho,rs, SNR_DB,rel_optDB);
       ee=norm(IM(:)'*x_t(:))/(IM(:)'*x_t(:));
       x_t_tv=ee*x_t;
       figure
       imshow(real(x_t_tv)*DeNSR/255)
       name_=[tv_name,'\Re_part_recons'];
       saveas(gcf,name_,'png')
       close

       figure
       imshow(imag(x_t_tv)*DeNSR/255)
       name_=[tv_name, '\Im_part_recons'];
       saveas(gcf,name_,'png')
       close
       
       
       
       noise_IM=x_t_-x_t_tv;
       figure
       imshow(real(noise_IM*DeNSR))
       name_=[tv_name,'\error_image'];
       saveas(gcf,name_,'png')
       close
       
       Rec_err(trial_SNR,trial_RHO,trial_RS)=-20*log(rel_opt)/log(10);  %rec_err(SNR_DB,RHO,rs)=err
       
       trial_RHO=trial_RHO+1; 
     end
     
     
     
     
     
     trial_RS=trial_RS+1;    
    end
        
        
   
trial_SNR=trial_SNR+1;             
end
%NSR_V(trial_RS,trial_SNR)

for i=1:len_RS
    KKK=zeros(len_RHO,l_d);
    
    
    for j=1:len_RHO
        for k=1:l_d
            
         KKK(j,k)=Rec_err(k,j,i);
        end
    end

[X,Y]=meshgrid(NSR_V(i,:),RHO);
figure
mesh(X,Y,KKK)
grid on
xlabel('SNR intensity')
ylabel('rho')
zlabel('SNR DB')
title_name=['rs=',num2str(RS(i))];
title(title_name)

loc=['/Users/Beaux/Desktop/phase retrieval/code/12/random_sample_chambolle/',title_name];
saveas(gcf,loc,'png')
    
    
end





