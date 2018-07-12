%%%%%% referee report about the paper : Solving large-scale general phase retrieval problems
%  via a sequence of convex relaxations

Bigname='/Users/Beaux/Desktop/phase retrieval'   ; % location of image

str1=[Bigname,'/image_lib/Cameraman.bmp'];   %% consider Im= Cameraman+i* barbera
str2=[Bigname,'/image_lib/Cameraman.bmp'];

os_rate=2;  %% oversampling ratio
TVswitch=0;    % TV   1==ON
            %      0==OFF
            


Noiseswitch = 1; % 1= add poisson noise to |Y|^2;        
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
MaxIter=40;
maxiter=8000;
Tol=1e-14;
toler=1e-15;  % tolerance in total variation denoising
tau= 1/8; % usually < 1/8;

SUBIM=6;


relative_DR_x=zeros(2,MaxIter);
DENSR=20:10:60;      %DENSR=15:20; % noise level used in the test
l_d=length(DENSR);

RHO=[0.5,0.8,1,1.3,10,50,100,600,1000];
rho = 4;
l_l= 0.2;
len_RHO=length(RHO);
NSR_V=zeros(1,length(DENSR));


% oversampling 
RS=0.7:0.1:1; % random sampling
len_RS=length(RS);
Rec_err=zeros(l_d,len_RHO,len_RS);   % reconstruction error rec_err(SNR_DB,RHO,rs)=err

[IM,Y1,Y2,SNR,mask]=proc_(str1,os_rate,num_mask,add_phase,...
                      sector_alpha,sector_beta,sector_mode,coer,SUBIM);
 if Noiseswitch == 0
     Y2=abs(Y1).^2;
 end
 

[Na,Nb]=size(IM);
norm_IM=norm(IM,'fro');
n=Na*Nb; %% dim of data vector
[Y_Na,Y_Nb]=size(Y1); %%% This Y_Nb= Nb*os_rate*'num_mask'
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

Bigname='/Users/Beaux/Desktop/phase retrieval/code/12/random_sample_chambolle';
mkdir(Bigname);


b=rand(Na,Nb); 
b = IM+0.0*rand(Na,Nb);
a=-b;

Ua=Nos_fft(a,os_rate,num_mask, mask);
Ub=Nos_fft(b,os_rate,num_mask, mask);
AA=Y2+conj(Ua).*Ua+conj(Ub).*Ua+conj(Ub).*Ub;
[l1,l2]=size(Ua);
    
X=[AA, conj(Ua+Ub); Ua+Ub,ones(size(Ua))];
post_n_norm=real(sum(AA(:)))+l1*l2;

rel_x=1;
rel_ADMM_x=zeros(MaxIter,1);



Y=zeros(size(X));

pre_n_norm=0;



copr=1;
while  rel_x > 1e-4    && copr < 2
     




     
     count=1;
     diff_n=1;
while diff_n > 1e-4  && count < 1000
    
    
    pre_n_norm=post_n_norm;
    
    Z=X+1/rho*Y;
    x=conj(Ub);
    Rex=real(x);
    Imx=imag(x);
    absx2=abs(x).^2;
    
    Re_Z1=real(Z(1:Y_Na,1:Y_Nb));
    Re_Z2=real(Z(1:Y_Na,1+Y_Nb:2*Y_Nb));
    Re_Z3=real(Z(1+Y_Na:2*Y_Na,1:Y_Nb));
    Im_Z2=imag(Z(1:Y_Na,1+Y_Nb:2*Y_Nb));
    Im_Z3=imag(Z(1+Y_Na:2*Y_Na,1:Y_Nb));
    ATuadmm_ucopr1=2*Rex.*(Re_Z1-Y2-absx2)+Re_Z2-Rex+Re_Z3-Rex;
    ATuadmm_ucopr2=2*Imx.*(Re_Z1-Y2-absx2)+Im_Z2-Imx+Im_Z3-Imx;
    ATuadmm_ucopr= [ATuadmm_ucopr1;ATuadmm_ucopr2];
    
    output=sqrt(2)*(ATuadmm_ucopr1.*Rex+ATuadmm_ucopr2.*Imx);
    Bx=1/2 *(ATuadmm_ucopr-sqrt(2)*[Rex.*(output./(ones(Y_Na,Y_Nb)+2*absx2)); Imx.*  (output./(ones(Y_Na,Y_Nb)+2*absx2)) ]);
    
    complex_Bx= Bx(1:Y_Na,1:Y_Nb)+1i*Bx(Y_Na+1:2*Y_Na,:);
    a_plus1=Nos_ifft(complex_Bx,os_rate,num_mask,mask);
    %%%lasso adjustment
    a_plus1=max(real(a_plus1)-l_l,0)+min(real(a_plus1)+l_l,0)+1i*(max(imag(a_plus1)-l_l,0)+min(imag(a_plus1)+l_l,0));
    
    % second stage in ADMM
    lambda=rho/2;
    
    Ua=Nos_fft(a_plus1,os_rate,num_mask, mask);
    Ub=Nos_fft(b,os_rate,num_mask, mask);
    AA=Y2+conj(Ua).*Ua+conj(Ub).*Ua+conj(Ub).*Ub;
    [l1,l2]=size(Ua);
    
    M=[AA, conj(Ua+Ub); Ua+Ub,ones(size(Ua))];
    
    
    C= M-1/rho*Y;
    C11=C(1:Y_Na,1:Y_Nb);
    C12=C(1:Y_Na,1+Y_Nb:2*Y_Nb);
    C21=C(1+Y_Na:2*Y_Na,1:Y_Nb);
    C22=C(1+Y_Na:2*Y_Na,1+Y_Nb:2*Y_Nb);
    lcc=length(C11(:));
    matrix_C= [spdiags(C11(:),0,lcc,lcc),spdiags(C12(:),0,lcc,lcc);spdiags(C21(:),0,lcc,lcc),spdiags(C22(:),0,lcc,lcc)];
    [U,D,V]=svd(full(matrix_C));
    D_X=max(0,D-1/lambda/2);
    X_plus1=U*D_X*V';
    %reshape XX
    reshape_Xplus1=[reshape(diag(X_plus1(1:Y_Na*Y_Nb,1:Y_Na*Y_Nb)),Y_Na,Y_Nb),reshape(diag(X_plus1(1:Y_Na*Y_Nb,1+Y_Na*Y_Nb:2*Y_Na*Y_Nb)),Y_Na,Y_Nb);...
        reshape(diag(X_plus1(1+Y_Na*Y_Nb:2*Y_Na*Y_Nb,1:Y_Na*Y_Nb)),Y_Na,Y_Nb),reshape(diag(X_plus1(1+Y_Na*Y_Nb:2*Y_Na*Y_Nb,1+Y_Na*Y_Nb:2*Y_Na*Y_Nb)),Y_Na,Y_Nb)];
    Y_plus1=Y+rho*(reshape_Xplus1-M);
    
    X=reshape_Xplus1;
    Y=Y_plus1;
    post_n_norm=real(sum(AA(:)))+l1*l2;
    diff_n=abs(post_n_norm-pre_n_norm);
    fprintf('count=%d, diff_n=%s\n',count, diff_n);
    count=count+1;
    
    
end 


b= - a_plus1;
ee=norm(IM(:)'*a_plus1(:))/(IM(:)'*a_plus1(:));
rel_x=norm(ee*a_plus1-IM,'fro')/norm(IM,'fro');
fprintf('COPR=%d, rel_x=%s\n', copr,rel_x)



rel_ADMM_x(copr)=rel_x;
copr=copr+1;
end
    
    
    
    
    

