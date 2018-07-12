% explore whether the 1-d image condition is necessary ? . not necessary.



Bigname='/Users/Beaux/Desktop/phase retrieval'   ; % location of image

str1=[Bigname,'/image_lib/Cameraman.bmp'];   %% consider Im= Cameraman+i* barbera
str2=[Bigname,'/image_lib/Cameraman.bmp'];



os_rate=2;  %% oversampling ratio


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
            
coer=3;    % regulating the SNR. coer higher, noise higher.
MaxIter=3000;
maxiter=2000;
Tol=1e-14;

SUBIM=150;


[IM,Y,Z,SNR,mask]=initial_1d_im(str1,str2,os_rate,num_mask,coer,SUBIM) ;
lambda_t= rand(size(Y));


if Noiseswitch==0
    Z=abs(Y).^2;
end

count_DR=1;
relative_DR_y =zeros(MaxIter,1);
rel_y=1;
gamma=1.2;


while rel_y > Tol && count_DR < MaxIter
    x_t=one_d_ifft(lambda_t,os_rate,num_mask,mask); 
    AAlambda=one_d_fft(x_t,os_rate,num_mask,mask);
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
   
            
    ee=norm(Y(:)'*AAlambda(:))/(Y(:)'*AAlambda(:));
    rel_y=norm(ee*AAlambda-Y,'fro')/norm(Y,'fro');
    relative_DR_y(count_DR)=rel_y;
    fprintf('count_DR=%d\n Coer=%f\n rel_y=%s\n', count_DR,DeNSR,rel_y);
    count_DR=count_DR+1;
               
end

figure
semilogy(1:count_DR-1, relative_DR_y(1:count_DR-1),'r');
xlabel('Iter')
ylabel('Err')
ee=norm(IM(:)'*x_t(:))/(IM(:)'*x_t(:));
a=real(ee*x_t); b = imag(ee*x_t);

imshow(b/255)




