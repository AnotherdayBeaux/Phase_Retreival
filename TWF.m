%%%%%% TWF  %%%%%%%






%%%% Compare modified-chenDR, DR_gau, DR_pois, TWF_pois




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
MaxIter=50;
Tol=1e-14;

SUBIM=50;

relative_DR_y=zeros(2,MaxIter);
DENSR=0.02:0.04:0.3; % noise level used in the test
NSR_V=zeros(1,length(DENSR));
end_error=zeros(4,length(DENSR));


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

DeNSR=0.05;

noise=normrnd(0,1,size(Y,1),size(Y,2));
NSR= norm(noise,'fro')/norm(IM,'fro');
B=abs(Y)+noise*DeNSR/NSR; % new abs(Y)
%%% verify NSR 
NSR= norm(B-abs(Y),'fro')/norm(Y,'fro');
fprintf('DeNSR=%f,\n NSR=%f\n',DeNSR,NSR);
%pause(2)
%NSR_V(trial)=NSR;

Z=B.^2;







z_0=IM.*rand(size(IM))+rand(size(IM))*15+rand(size(IM))*1i*40;


a_up=6; a_low=0.05; a_h=7;
            [M1,M2]=size(Y);
        norm_a=zeros(size(Y));
if aaaa==0
        for i=1:M1
            for j=1:M2
            Ind=zeros(size(Y));
            Ind(i,j)=1;
            aij=Nos_ifft(Ind,os_rate,num_mask,mask);
            norm_a(i,j)=norm(aij,'fro');

            end
        end
      sqrt_n= sqrt(size(z_0,1)*size(z_0,2));
end      
    %%%% TWF 
    z_t=z_0;
    count_TWF=1;
while count_TWF <30000



        fft_z_t=Nos_fft(z_t,os_rate,num_mask,mask);
        abs_fft_z_t=abs(fft_z_t);
        eps1=sqrt_n*(abs_fft_z_t./norm_a)/norm(z_t,'fro');

        Ind_1a=zeros(size(Y));
        Ind_1b=zeros(size(Y));
        
        Ind_1a(find(a_low<eps1))=1;
        Ind_1b(find(a_up>eps1))=1;
        Ind_1=Ind_1a.*Ind_1b;


        %eps2

        left_eps2=abs(Z-abs_fft_z_t.^2);
        K_t=sum(sum(left_eps2))/size(B,1)/size(B,2);

        Ind_2=zeros(size(Y));
        Ind_2(find(left_eps2<K_t*a_h*eps1))=1;


        multi_d_orig=(Z-abs_fft_z_t.^2)./conj(fft_z_t);
        multi_d=multi_d_orig.*Ind_1.*Ind_2;
        grad_t=2/size(Y,1)/size(Y,2)*Nos_ifft(multi_d,os_rate,num_mask,mask);
        grad=2/size(Y,1)/size(Y,2)*Nos_ifft(multi_d_orig,os_rate,num_mask,mask);
        mu=1; % step length



        l_z_pre=1;
        l_z_post=0;
        v_incre=abs(grad.*conj(grad));% modify grad_t
        incre=sum(sum(v_incre));
        
        while l_z_pre > l_z_post
    
            A_z=Nos_fft(z_t+grad*mu,os_rate,num_mask,mask); % modify grad_t
            ll_z_post=2*(Z).*log(abs(A_z))-abs(A_z).^2;
            l_z_post=sum(sum(ll_z_post));
            A_zz=Nos_fft(z_t,os_rate,num_mask,mask);
            ll_z_pre=2*(Z).*log(abs(A_zz))-abs(A_zz).^2;
            l_z_pre=sum(sum(ll_z_pre))+mu*incre;
            
    
    
            mu=0.2*mu; 
            
        end
        mu=mu*5;
        z_t=z_t+grad*mu;% modify grad_t
        
        
        ee=norm(IM(:)'*z_t(:))/(IM(:)'*z_t(:));
        rel_x=norm(ee*z_t-IM,'fro')/norm(IM,'fro');
        
        fprintf('count=%d\n rel_x=%f\n',count_TWF,rel_x)
        
        count_TWF=count_TWF+1;
        
        
        
        
        
                     
end



