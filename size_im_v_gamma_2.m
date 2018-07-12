%%%%%% implement DR method to poisson likelihood function


str='/Users/Beaux/Desktop/phase retrieval/image_lib/Cameraman.bmp';

os_rate=2;  %% oversampling ratio
TVswitch=0;    % TV   1==ON
            %      0==OFF
tv=0.1; % tv * ||grad_image||_1
Noiseswitch = 0; % 1= add poisson noise to |Y|^2;        
num_mask=3; % 0  0 mask / Id mask
            % 1  1 mask case
            % 2  2 mask case 
            % 3  1 and 1/2 mask case
            
SUB_IM=56:20:256;
Gamma=0.4:0.2:30;
compare_iter=zeros(length(56:20:256),length(Gamma));
MaxIter=10000;
Tol=1e-15;
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







num_image=1;
for SUBIM=SUB_IM
    
    
    
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

trial=1;
%%% 3 way initialization
%1
%lambda_0=Nos_fft(Nos_ifft(sqrt(Z),os_rate,num_mask,mask),os_rate,num_mask,mask);
%2
lambda_0=rand(size(Y)).*exp(2i*pi*rand(size(Y)));
%3
%lambda_0=Y+rand(size(Y))*1e-1;

for gamma=Gamma
    
    
lambda_t=lambda_0;
y_tplus1=zeros(size(lambda_t));
count=1;  
residual=zeros(MaxIter,1);

resi=1;

while resi > Tol && count < MaxIter
    
    %%%%%% y_t+1 = argmin 1/2||y-A'A\lambda_t||^2+1/2gamma ||y-lambda_t||^2 %%%%%%%%
    x_t=Nos_ifft(lambda_t,os_rate,num_mask,mask);
    AAlambda=Nos_fft(x_t,os_rate,num_mask,mask);
    %y_tplus1=(gamma*AAlambda+lambda_t)/(gamma+1);
    y_tplus1=AAlambda;
    
    %%%%% z_t+1 = argmin |z|^2-b log(|z|^2)+1/2gamma*||2y_t+1-lambda_t-z||^2
    ylambda_k=2*y_tplus1-lambda_t;
    Q=abs(ylambda_k);
    rot_z=(Q/gamma+ sqrt(Q.^2/gamma^2+8*(2+1/gamma)*Z))/(4+2/gamma);
    ang_z_real=real(ylambda_k)./Q;
    ang_z_imag=imag(ylambda_k)./Q;

    z_tplus1=rot_z.*ang_z_real+1i*rot_z.*ang_z_imag;
    
    %%%% lambda_t+1=lambda_t+z_t+1-y_t+1;
    
    lambda_t=lambda_t+z_tplus1-y_tplus1;
    
    resi=norm(abs(y_tplus1)-sqrt(Z),'fro')/norm(sqrt(Z),'fro');
    residual(count)=resi;
    %fprintf('count=%d\n resi=%f\n', count,resi)
    count=count+1;
    
end
%Name=['/Users/Beaux/Desktop/phase retrieval/code/3/',...
%    'add_phase=',num2str(add_phase),'_gamma=',num2str(gamma*100),...
%    'Noiseswitch=',num2str(Noiseswitch)];

%mkdir(Name)



%figure
%semilogy(1:count-1,residual(1:count-1),'g')
%title('Relative error')
%xlabel('Iter')
%ylabel('|||A^{*}Alambda|-b||_2/||b||_2')
%name=[Name,'/relative_error'];
%saveas(gcf,name,'png')
%close


compare_iter(num_image,trial)=count;
trial=trial+1
end
num_image=num_image+1

end
figure
[X,Y]=meshgrid(SUB_IM,Gamma);
surf(X,Y,compare_iter)
ylabel('gamma')
xlabel('size of image n-by-n')
name=['/Users/Beaux/Desktop/phase retrieval/code/3/','opt_gamma'];
saveas(gcf,name,'png')
close

plot(Gamma,compare_iter(10,:),'g')
hold on
plot(Gamma,compare_iter(9,:),'r')
hold on
plot(Gamma,compare_iter(8,:),'b')
hold on
plot(Gamma,compare_iter(7,:),'k')
hold on
plot(Gamma,compare_iter(6,:),'g:')
hold on
plot(Gamma,compare_iter(5,:),'r:')
hold on
plot(Gamma,compare_iter(4,:),'b:')
hold on
xlabel('Gamma 0.4:0.2:30')
ylabel('iter taken to 1e-15')
legend('236-by-236','216-by-216','196-by-196','176-by-176',...
    '156-by-156','136-by-136','116-by-116')

Opt=zeros(10,1);
for i = 1:10
    
    aaa=find(compare_iter(i,:)==min(compare_iter(i,:)));
    opt_gamma=Gamma(aaa(1));
    Opt(i)=opt_gamma;
 
end

plot(SUB_IM(1:10),Opt,'gs-')
xlabel('size of IM n-by-n')
ylabel('optimal gamma')
axis([40 256 0 30])







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
