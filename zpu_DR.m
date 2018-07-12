%%%   |z|^2-b log(|z|^2)+alpha||p||_1 + DR 
%%%   sub to p=grad u; y= A^*u 
%%%





% pois likelihood 

% gaus likelihood 



Bigname='/Users/Beaux/Desktop/phase retrieval'   ; % location of image


str1=[Bigname,'/image_lib/Cameraman.bmp'];   %% consider Im= Cameraman+i* barbara
str2=[Bigname,'/image_lib/barbara.png'];

os_rate=2;  %% oversampling ratio = 4
rs=1; % rs=downsampling ratio.
% the total sample points collected = rs*os_rate^2

Noiseswitch = 0; % 1= add poisson noise to |Y|^2;        
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

SUBIM=156;
type_im=0;  % 0 real
            % 1 str1+i*str1
            % 2 str1+i*str2
            
type_f=1;   % 1 pois
            % 0 gaus
            
type_l1=0;  % 1 lengthwise
            % 0 component wise



relative_DR_x=zeros(1,MaxIter);
DENSR=20:10:60;      %DENSR=15:20; % noise level used in the test
l_d=length(DENSR);

RHO2=[5,50,200,1000,5000,10000,50000];
RS=0.1:0.2:1;
%RHO2=2000;
RS=1;
RHO2=5000;
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


 





Bigname='/Users/Beaux/Desktop/phase retrieval/code/15/DR_test';
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
    
    q_0=grad(IM);
    re_q_0=real(q_0);
    im_q_0=imag(q_0);
    re_q_0(re_q_0>0)=re_q_0(re_q_0>0)+thresh;
    re_q_0(re_q_0<0)=re_q_0(re_q_0<0)-thresh;
    
    im_q_0(im_q_0>0)=im_q_0(im_q_0>0)+thresh;
    im_q_0(im_q_0<0)=im_q_0(im_q_0<0)-thresh;
    q_0=im_q_0*1i+re_q_0;
    
    
    
    
    s_0=Nos_fft(IM,os_rate,num_mask,mask);
    q_0=rand([size(IM),2]);
    s_0=rand(size(Y));
    
    
    
   
    
    
    for rho2=RHO2
        %DR_Pois
        q_k=q_0;
        s_k=s_0;

       
        
        
        count_DR=1; 
        rel_up=1e5;
        rel_u=1e4;
        rho_lap=spdiags(rho2/rho1*ones(Na*Nb,1),0,Na*Nb,Na*Nb)-lap;
        
        
        while (rel_up > rel_u && count_DR < MaxIter) || count_DR<100
            J_G_q=J_G(rho1,alpha,q_k,type_l1);
            J_F_s=J_F(Z,rho2,s_k,RanS,type_f);
            updated_s=s_k-2*J_F_s;
            inter=div(q_k-2*J_G_q)-rho2/rho1*Nos_ifft(updated_s,os_rate,num_mask,mask);
            
         
            b=rho_lap\inter(:);
            u_kplus1=reshape(b,Na,Nb);
            
            q_kplus1=grad(u_kplus1)+q_k-J_G_q;
            s_kplus1=Nos_fft(u_kplus1,os_rate,num_mask,mask)+s_k-J_F_s;
            
            %diff_q_k=norm(q_k-q_kplus1,'fro')/norm(q_k,'fro');
            diff_s_k=norm(s_k-s_kplus1,'fro')/norm(s_k,'fro');
            q_k=q_kplus1;
            s_k=s_kplus1;
            u_k_=Nos_ifft(s_k,os_rate,num_mask,mask);
            
            %ee=norm(IM(:)'*u_kplus1(:))/(IM(:)'*u_kplus1(:));
            ee=norm(IM(:)'*u_k_(:))/(IM(:)'*u_k_(:));
            rel_up=rel_u;
            %rel_u=norm(ee*u_kplus1-IM,'fro')/norm(IM,'fro');
            rel_u=norm(ee*u_k_-IM,'fro')/norm(IM,'fro');
            relative_DR_x(count_DR)=rel_u;
            %fprintf('count_DR=%d\n DB_int=%f\n rel_u=%s\n', count_DR,SNR_DB,rel_u);
            fprintf('count_DR=%d\n d_s_k=%f\n rel_u=%s\n', count_DR,diff_s_k,rel_u);
            count_DR=count_DR+1;
            
        end
        Rec_err(trial_SNR,trial_RHO2,trial_RS)=-10* log(rel_u^2)/log(10);
        
        
        ee=norm(IM(:)'*u_k_(:))/(IM(:)'*u_k_(:));
        u_k_=ee*u_k_;
        figure
        imshow(real(u_k_)*DeNSR/255)
        figure
        imshow(imag(u_k_)*DeNSR/255)
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