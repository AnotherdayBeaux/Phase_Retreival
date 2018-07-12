
% phitography






%%%   noiseless case   DR or ADMM 

%%%   noise     case   
%%%





% pois likelihood 

% gaus likelihood 





Bigname='/Users/Beaux/Desktop/phase retrieval'   ; % location of image


str1=[Bigname,'/image_lib/Cameraman.bmp'];   %% consider Im= Cameraman+i* barbara
str2=[Bigname,'/image_lib/barbara.png'];

os_rate=2;  %% oversampling ratio = 4
rs=1; % rs=downsampling ratio.
% the total sample points collected = rs*os_rate^2

Noiseswitch = 1; % 1= add poisson noise to |Y|^2;        
num_mask=4; % 0  0 mask / Id mask
            % 1  1 mask case
            % 2  2 mask case 
            % 3  1 and 1/2 mask case
            % 4  1 id, 2 random ones.
            
                         
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

add_phase=0;% 0 real image
            % 1 add phase consistent with sector constrain
            
coer=3;    % regulating the SNR. coer higher, noise higher.
MaxIter=20000;
%maxiter=2000;
Tol=1e-10;


% tv denoising parameter. 
%toler=1e-12;  % tolerance in total variation denoising
%tau= 1/10;    % usually < 1/8;
%tv=0.1;

SUBIM=256;
type_im=2;  % 0 real
            % 1 str1+i*str1
            % 2 str1+i*str2


relative_DR_x=zeros(1,MaxIter);
DENSR=20:10:60;      %DENSR=15:20; % noise level used in the test
DENSR=60;
l_d=length(DENSR);

RHO2=[5,50,200,1000,5000,10000,50000];
RHO2=5000;
RS=0.1:0.2:1;
RS=0.2;
%RHO2=2000;

len_RS=length(RS);
len_RHO2=length(RHO2);
NSR_V=zeros(len_RS,l_d);
Rec_err=zeros(l_d,len_RHO2,len_RS);   % reconstruction error

[IM,Y,Z,SNR_DB,mask]=proc_(str1,os_rate,num_mask,add_phase,...
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
%w1w2Matrix=16*((sin(pi*(0:Na-1)/Na)).^2)'*(sin(pi*(0:Nb-1)/Nb)).^2;
B=[ones(Na*Nb,1),ones(Na*Nb,1),-4*ones(Na*Nb,1),ones(Na*Nb,1),ones(Na*Nb,1)];

lap=spdiags(B,[-Na,-1,0,1,Na],Na*Nb,Na*Nb);




% initialize  lambda1_k lambda2_k lambda3_
%             rho1    rho2    rho3
%             y_k     p_k     z_k
%             u_k
%
%%%%%%%%%%%%%%%%%%%%%%
rho1=2;    %%%%%%%
alpha=200; %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
thresh=alpha/rho1; % update p
lambda2_0=rand(size(Y)).*exp(1i*rand(size(Y))*pi)*0;
lambda1_0=rand([size(IM),2]).*exp(1i*rand([size(IM),2])*pi)*0;


 





Bigname='/Users/Beaux/Desktop/phase retrieval/code/14/random_sample';
mkdir(Bigname);
trial_SNR=1;

%%% add poisson random noise
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
    NSR= norm(B-abs(Y).*RanS,'fro')/norm(Y.*RanS,'fro');
    SNR_DB= -10* log(NSR^2)/log(10);
    fprintf('DeNSR=%f,\n NSR=%f\n SNR_DB=%f\n',DeNSR,NSR,SNR_DB);

    NSR_V(trial_RS,trial_SNR)=SNR_DB;
    %%%% end verify   
 
    %%%%% FOR GAMMA
    name=[Bigname,'/DeNSR=',num2str(DeNSR),',rs=',num2str(rs)];   %% modify 'DeNSR' to 'NSR' 
    % when noiseswitch is changed to 1.
    mkdir(name)
    
    trial_RHO2=1;
    for rho2=RHO2
        %DR_Pois
        lambda1_k=lambda1_0;
        lambda2_k=lambda2_0;

        
        u_k=Nos_ifft(RanS.*B,os_rate,num_mask,mask);
        z_k=RanS.*B;
        p_k=grad(u_k);
        
        
        count_DR=1; 
        rel_up=1e5;
        rel_u=1e4;
        
        rho_lap=spdiags(rho2/rho1*ones(Na*Nb,1),0,Na*Nb,Na*Nb)-lap;
        
        
        while (rel_up > rel_u && count_DR < MaxIter) || count_DR<100
            
            %%%% 
            % update u
            
             combo_z_lambda2=z_k+lambda2_k/rho2;
             combo_p_lambda1=p_k+lambda1_k/rho1;
             Azlambda2=Nos_ifft(combo_z_lambda2,os_rate,num_mask,mask);
             divplambda1=div(combo_p_lambda1);
             Right=rho2/rho1*Azlambda2-divplambda1;
             b=rho_lap\Right(:);
             u_kplus1=reshape(b,Na,Nb);
             
             
             
             %fft2_co=fft2(rho2/rho1*Azlambda2-divplambda1);
             %u_kplus1=ifft2(fft2_co./(rho2/rho1+w1w2Matrix));
            
            
            
            % update p
            grad_u=grad(u_kplus1);
            inter=grad_u-lambda1_k/rho1;
            re_inter=real(inter);
            im_inter=imag(inter);
            
            re_inter(re_inter>0)=1;
            re_inter(re_inter<0)=-1;
            
            im_inter(im_inter>0)=1;
            im_inter(im_inter<0)=-1;  
            
            real_inter=real(inter);
            imag_inter=imag(inter);
            
            p_kplus1=max(0,abs(real_inter)-thresh).*re_inter+1i*max(0,abs(imag_inter)-thresh).*im_inter;
            
            lambda1_kplus1=lambda1_k+rho1*(p_kplus1-grad_u);   % update lambda1

            % update z
            Au=Nos_fft(u_kplus1,os_rate,num_mask,mask);
            
            combo_Au_lambda2=Au-lambda2_k/rho2;
            Q=abs(combo_Au_lambda2);
            % POISSON
            rot_z=rho2/(4+2*rho2)*(Q +sqrt(Q.^2+8*Z*(2+rho2)/rho2^2));
            
            % GAUSSIAN
            %rot_z= (sqrt(Z)+rho2*Q)/(1+rho2);
            
            z_kplus1=rot_z.*combo_Au_lambda2./Q;
            
            z_kplus1(isnan(z_kplus1))=0;
            z_kplus1=z_kplus1.*RanS;
            z_kplus1=z_kplus1.*RanS+abs(RanS-1).*combo_Au_lambda2;
            
            lambda2_kplus1=lambda2_k+rho2*(z_kplus1-Au);    % update lambda2      
   
            
            
            ee=norm(IM(:)'*u_kplus1(:))/(IM(:)'*u_kplus1(:));
            rel_up=rel_u;
            rel_u=norm(ee*u_kplus1-IM,'fro')/norm(IM,'fro');
            relative_DR_x(count_DR)=rel_u;
            fprintf('count_DR=%d\n DB_int=%f\n rel_u=%s\n', count_DR,SNR_DB,rel_u);
            
            count_DR=count_DR+1;
            
            lambda1_k=lambda1_kplus1;
            lambda2_k=lambda2_kplus1;
        
            z_k=z_kplus1;
            u_k=u_kplus1;
            p_k=p_kplus1;
        end
        Rec_err(trial_SNR,trial_RHO2,trial_RS)=-10* log(rel_u^2)/log(10);
        
        
        ee=norm(IM(:)'*u_k(:))/(IM(:)'*u_k(:));
        u_k_=ee*u_k;
        figure
        imshow(real(u_k_)*DeNSR/255)
        imshow(imag(x_t_)*DeNSR/255)
        name_=[name,'/SNR_DB=',num2str(SNR_DB),'_rs=',num2str(rs),'_rho2=',num2str(rho2)];
        mkdir(name_)
        
        name_im=[name_,'/recon_im'];
        %saveas(gcf,name_im,'png')
        %close
        figure
        semilogy(1:count_DR-1, relative_DR_x(1:count_DR-1),'g')
        xlabel('interation')
        ylabel('rel_u')
        Title=['SNR DB=',num2str(SNR_DB),',rs=',num2str(rs),',rho2=',num2str(rho2)];
        title(Title)
        
        name_plot=[name_,'/rel_u_plot'];        
        %saveas(gcf,name_plot,'png')
        %close
        
    trial_RHO2=trial_RHO2+1;    
    end
    trial_RS=trial_RS+1;
    
end



trial_SNR=trial_SNR+1;

end