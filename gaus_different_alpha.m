%%%% consider the gaussian likelihood function in the case when alpha = 1/2




Bigname='/Users/Beaux/Desktop/phase retrieval'   ; % location of image

str1=[Bigname,'/image_lib/Cameraman.bmp'];   %% consider Im= Cameraman+i* barbera
str2=[Bigname,'/image_lib/Cameraman.bmp'];

os_rate=2;  %% oversampling ratio


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
            
coer=3;    % regulating the SNR. coer higher, noise higher.
MaxIter=3000;
maxiter=2000;
Tol=1e-14;
eta=-1/2+sqrt(3)/2 * 1i;
SUBIM=100;

relative_DR_y=zeros(1,MaxIter);
DENSR=1:4:20;      %DENSR=15:20; % noise level used in the test
l_d=length(DENSR);

RHO=[0.05,0.1,10,20];
Gamma=1./RHO;       %Gamma=[1, 2.5,4];
l_g=length(Gamma);
NSR_V=zeros(1,length(DENSR));
end_error=zeros(2*l_g,length(DENSR));

Rec_err=zeros(l_d,l_g);   % reconstruction error rec_err(nsr,gamma,tv)=err

[IM,Y,Z,SNR,mask]=proc_TV(str1,str2,os_rate,num_mask,add_phase,...
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

Bigname='/Users/Beaux/Desktop/phase retrieval/code/17/different alpha';
mkdir(Bigname);
trial=1;


for DeNSR=DENSR
    [IM,Y,Z,SNR,mask]=proc_TV(str1,str2,os_rate,num_mask,add_phase,...
                      sector_alpha,sector_beta,sector_mode,DeNSR,SUBIM);
if Noiseswitch == 0
     Z=abs(Y).^2;
end
B=sqrt(abs(Z));  %%% noise case
              
%%% verify NSR 
NSR= norm(B-abs(Y),'fro')/norm(Y,'fro');
fprintf('DeNSR=%f,\n NSR=%f\n',DeNSR,NSR);
NSR_V(trial)=NSR;
%%%% end verify   

    %lambda_0=Y+rand(size(Y));
    lambda_0=rand(size(Y))+1i*rand(size(Y));
    gamma_count=1;
    
    
    
    %%%%% FOR GAMMA
    name=[Bigname,'/NSR=',num2str(NSR)];   %% modify 'DeNSR' to 'NSR' 
    % when noiseswitch is changed to 1.
    mkdir(name)
    
    
    for gamma=Gamma
    %DR_Pois
    lambda_t=lambda_0;
    count_DR=1;   
    rel_y=1;
    while rel_y >Tol  && count_DR < MaxIter
    
            x_t=Nos_ifft(lambda_t,os_rate,num_mask,mask); 
            AAlambda=Nos_fft(x_t,os_rate,num_mask,mask);
            y_tplus1=AAlambda;

            ylambda_k=2*y_tplus1-lambda_t;
            Q=abs(ylambda_k);
            % POISSON
            %rot_z=(Q/gamma+sqrt(Q.^2/gamma^2+8*(2+1/gamma)*(Z)))/(4+2/gamma);
            
            % GAUSSIAN
            f = @(z) 1/2*(sqrt(z)-sqrt(B)).^2+1/gamma/2*(z-Q).^2;
            q = -gamma/2* sqrt(B);
            p= +gamma/2 - Q;
            
            delta= (q.^2)/4+ (p.^3)/27;
            ind_neg=zeros(size(delta));   % ind_neg pick out the case when delta <0.
            ind_neg(delta<0)=1;
            ind_pos=zeros(size(delta));   % 
            ind_pos(delta>0)=1;
            
             a=(-q/2+sqrt(q.^2/4+p.^3/27)).^(1/3);
             b_pre=-q/2-sqrt(q.^2/4+p.^3/27);
             b=abs(-q/2-sqrt(q.^2/4+p.^3/27)).^(1/3);
             b=b.*(2*(b_pre>0)-1);
             conj_a= conj(a);
             root1=(a+conj_a).*ind_neg;
             root2=(eta*a+conj(eta)*conj_a).*ind_neg;
             root3=(conj(eta)*a+eta*conj_a).*ind_neg;
             root1=max(root1,0);
             root2=max(root2,0);
             root3=max(root3,0);
             
             f_root1=f(root1);
             f_root2=f(root2);
             f_root3=f(root3);
             
            
             minfroot=min(f_root1,f_root2);
             minfroot=min(minfroot,f_root3);
             
             
             
            rot_z= ((a+conj_a).*ind_neg+ (a+b).*ind_pos ).^2;
            
            ang_z_real=real(ylambda_k)./Q;
            ang_z_imag=imag(ylambda_k)./Q;
            z_tplus1=rot_z.*ang_z_real+1i*rot_z.*ang_z_imag;
    
            lambda_t=lambda_t+z_tplus1-y_tplus1;
   
            
            ee=norm(Y(:)'*AAlambda(:))/(Y(:)'*AAlambda(:));
            rel_y=norm(ee*AAlambda-Y,'fro')/norm(Y,'fro');
            relative_DR_y(count_DR)=rel_y;
            fprintf('count_DR=%d\n NSR=%f\n rel_y=%s\n', count_DR,NSR,rel_y);
            count_DR=count_DR+1;

    end
    ee=norm(IM(:)'*x_t(:))/(IM(:)'*x_t(:));
    x_t_=ee*x_t;
    
    Title=['gamma=',num2str(gamma),'; NSR=',num2str(NSR)];
    figure
    %imshow(real(x_t_)*DeNSR/255)
    imshow(real(x_t_)/255)
    name_imag=[name,'/',Title,',image'];
    saveas(gcf,name_imag,'png')
    close
    
    figure
    semilogy(1:count_DR-1,relative_DR_y(1:count_DR-1),'r');
    xlabel('iter')
    ylabel('rel err')
    Title=['gamma=',num2str(gamma),'; NSR=',num2str(NSR)];
    title(Title)
    name_iter=[name,'/',Title];
    saveas(gcf,name_iter,'png')
    close

    %%%%%%%%%%%%%
    %%%%%%%%%%%%%
       

     
     
     
     
     Rec_err(trial,gamma_count) =rel_y ;  
     gamma_count=gamma_count+1;    
    end
    
 trial =trial+1;   
end
  figure  
for i=1:4
    
    NSR=NSR_V(i);
    
    plot(Gamma,Rec_err(i,:));
    hold on
    
    
end
NSR1=['NSR=',num2str(NSR_V(1))];
NSR2=['NSR=',num2str(NSR_V(2))];
NSR3=['NSR=',num2str(NSR_V(3))];
NSR4=['NSR=',num2str(NSR_V(4))];
legend(NSR1,NSR2,NSR3,NSR4);
xlabel('gamma')
ylabel('rec err')

[X,Y]=meshgrid(NSR_V,Gamma);
figure
mesh(X,Y,Rec_err')
xlabel('NSR')
ylabel('Gamma')
zlabel('Rec err')