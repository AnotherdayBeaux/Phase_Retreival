%%%% magnitude retreival 


%%%% when m <2n+1 of course no uniqueness.
%%%% when m=2n+1 ?? 

% still use DR method. the update becomes 

% y update: y_t=A^*A lambda_t
% z update: z_t= abs(2y_t-lambda_t)\odot phase(2y_t-lambda_t) 
% lambda update lambda_t+1=lambda_t+ z_t-y_t.

Bigname='/Users/Beaux/Desktop/phase retrieval'   ; % location of image

str=[Bigname,'/image_lib/Cameraman.bmp'];

os_rate=2;  %% oversampling ratio
rs =0.7;     % random sampling ratio.
TVswitch=0;    % TV   1==ON
            %      0==OFF
tv=0.1; % tv * ||grad_image||_1
Noiseswitch = 1; % 1= add poisson noise to |Y|^2;        
num_mask=2; % 0  0 mask / Id mask
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
MaxIter=2000;
Tol=1e-14;

SUBIM=156;

relative_DR_y=zeros(1,MaxIter);
%Gamma=[0.5,1,1.5,2.5,4];
Gamma=1;
l_g=length(Gamma);



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

% Add random sampling 
RANS=rand(1,Y_Na*Y_Nb);
pos_1=find(RANS<rs);
RANS=zeros(1,Y_Na*Y_Nb);
RANS(pos_1)=1;
RanS=reshape(RANS,Y_Na,Y_Nb);
trs=Y_Na*Y_Nb/sum(RanS(:));
phase=Y./abs(Y);



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

lambda_0=rand(size(Y));
    gamma_count=1;
    for gamma=Gamma 
   
    lambda_t=lambda_0;
    count_DR_gau=1;   
    resi=1;
    while resi > Tol && count_DR_gau < MaxIter
                
        
            %lambda_t = RanS.* lambda_t; % rs 
            x_t=Nos_ifft(lambda_t,os_rate,num_mask,mask);
            AAlambda=Nos_fft(x_t,os_rate,num_mask,mask);
            y_tplus1=AAlambda;
            %y_tplus1=RanS.*y_tplus1;   % rs
      
            ylambda_k=2*y_tplus1-lambda_t;
            Q=abs(ylambda_k);
            phase_ylk=ylambda_k./Q;
            hori_ylk=phase_ylk.*conj(phase);
            hori_ylk(isnan(hori_ylk))=0;
            re_ylk=real(hori_ylk);
            phase_modu=ones(size(re_ylk));
            phase_modu(re_ylk<0)=-1;
            phase_modu(re_ylk>0)=1;
            
            z_tplus1=Q.*phase.*phase_modu.*RanS+Q.*(1-RanS).*phase_ylk;
    
            lambda_t=lambda_t+z_tplus1-y_tplus1;
            
            x_t=Nos_ifft(ylambda_k,os_rate,num_mask,mask);
            nor_x_t=x_t*norm(IM,'fro')/norm(x_t,'fro');
            rel_y=norm(nor_x_t-IM,'fro')/norm(IM,'fro');
            relative_DR_y(1,count_DR_gau)=rel_y;
            fprintf('count=%d\n rel_y=%f\n', count_DR_gau,rel_y);
            count_DR_gau=count_DR_gau+1;
    end
    hold on
      loglog(1:count_DR_gau-1,relative_DR_y(1:count_DR_gau-1),'r')
      xlabel('iter')
      ylabel('rel im')
      
      gamma_count=gamma_count+1;
        
    end