%%%%%% implement DR method to poisson likelihood function
%%%%%% Assume the gamma in the 1st, 2nd updates are different.

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

SUBIM=256;
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
%%% 3 way initialization
%1
%lambda_0=Nos_fft(Nos_ifft(sqrt(Z),os_rate,num_mask,mask),os_rate,num_mask,mask);
%2
lambda_0=rand(size(Y)).*exp(2i*pi*rand(size(Y)));
%3
%lambda_0=Y+rand(size(Y))*1e-1;
MaxIter=10000;
Tol=1e-12;
trial_ga1=0;
Gamma_1st=5:100;
Gamma_2nd=3:50;
compare_iter=zeros(length(Gamma_1st),length(Gamma_2nd));
compare_lambda_2=zeros(length(Gamma_1st),length(Gamma_2nd));
Bigname='/Users/Beaux/Desktop/phase retrieval/code/2/fig_nonoise_Cameraman_compare_ga1_ga2';
mkdir(Bigname)

for gamma_1st = Gamma_1st
    %%%%%%%%%%%
    trial_ga1=trial_ga1+1; % count in gamma_1st
    trial_ga2=0;
    %%%%%%%%%%%
    for gamma_2nd= Gamma_2nd
        %%%%%%%%
        trial_ga2=trial_ga2+1;   % count in gamma_2nd
        %%%%%%%%  
        lambda_t=lambda_0;
        y_tplus1=zeros(size(lambda_t));
        count=1;  
        residual=zeros(MaxIter,1);
        resi1=zeros(MaxIter,1);
        %dev=zeros(MaxIter,1);

        resi=1;

        while resi > Tol && count < MaxIter
    
            %%%%%% y_t+1 = argmin 1/2||y-A'A\lambda_t||^2+1/2gamma ||y-lambda_t||^2 %%%%%%%%
            x_t=Nos_ifft(lambda_t,os_rate,num_mask,mask);
            AAlambda=Nos_fft(x_t,os_rate,num_mask,mask);
            %Resi=norm(y_tplus1,'fro')^2+norm(AAlambda,'fro')^2-2*abs(sum(sum(y_tplus1.*conj(AAlambda))));
            %resi1(count)=sqrt(norm(Resi))/norm(y_tplus1,'fro');
            y_tplus1=(gamma_1st*AAlambda+lambda_t)/(gamma_1st+1);
            %y_tplus1=AAlambda;
      
    
            %%%%% z_t+1 = argmin 1/2 |||z|-b||^2 +1/2gamma*||2y_t+1-lambda_t-z||^2
            ylambda_k=2*y_tplus1-lambda_t;
            Q=abs(ylambda_k);
            %rot_z=(Q/gamma_2nd+ sqrt(Q.^2/gamma_2nd^2+8*(2+1/gamma_2nd)*Z))/(4+2/gamma_2nd);
            rot_z=(Q/gamma_2nd+sqrt(Z))/(1+1/gamma_2nd);
            ang_z_real=real(ylambda_k)./Q;
            ang_z_imag=imag(ylambda_k)./Q;

            z_tplus1=rot_z.*ang_z_real+1i*rot_z.*ang_z_imag;
    
            %%%% lambda_t+1=lambda_t+z_t+1-y_t+1;
    
            lambda_t=lambda_t+z_tplus1-y_tplus1;
    
            resi=norm(abs(y_tplus1)-sqrt(Z),'fro')/norm(sqrt(Z),'fro');
            residual(count)=resi;
    
            %kkkkk=norm(x_t,'fro')^2+norm_IM^2-2* abs(sum(sum(IM.*conj(x_t))));
            %dev(count)=sqrt(norm(kkkkk))/norm_IM;
            count=count+1;
    
        end
        Resi1=norm(y_tplus1,'fro')^2+norm(AAlambda,'fro')^2-...
            2*abs(sum(sum(y_tplus1.*conj(AAlambda))))/norm(y_tplus1,'fro');
        fprintf('count=%d\n gamma_1st=%f\n gamma_2nd=%f\n', count,gamma_1st,gamma_2nd)
        Name=[Bigname,...
            '/add_phase=',num2str(add_phase),'_ga1st=',...
            num2str(gamma_1st),'_ga2nd=',num2str(gamma_2nd),...
            'Noiseswitch=',num2str(Noiseswitch)];

        mkdir(Name)
        conf=[Name,'/configration.txt'];
        fid=fopen(conf,'wt');
        fprintf(fid, 'gamma1st= %f\n gamma2nd=%f\n Noise=%d\n number of iter=%d\n Tol=%f\n n end_resi=%f\n end_lambda_2=%f\n'...
                 , gamma_1st,gamma_2nd, Noiseswitch, (count-1),Tol,resi,Resi1);
        fclose(fid);



    figure
    semilogy(1:count-1,residual(1:count-1),'g')
    title('semily Relative error')
    xlabel('Iter')
    ylabel('|||A^{*}Ax|-b||_2/||b||_2')
    name=[Name,'/relative_error'];
    saveas(gcf,name,'png')
    close
   



%figure
%plot(count-1000:count-1,residual(count-1000:count-1),'g')
%title('last 1000 relative error (Magnified)')
%xlabel('Iter')
%ylabel('|||A^Ax|-b||_2/||b||_2')
%name=[Name,'/relative_error_last50'];
%saveas(gcf,name,'png')
%close

%%%%% 2nd choice of stopping criteria
    compare_iter(trial_ga1,trial_ga2)=count;
    compare_lambda_2(trial_ga1,trial_ga2)=Resi1;
    end
end
figure
[X,Y]=mesh(Gamma_1st,Gamma_2nd);
surf(X,Y,compare_iter)
xlabel('gamma_1st=')
ylabel('gamma_2nd=')
zlabel('iter to reach 1^{-15}')
name=[Bigname,'opt_gamma'];
saveas(gcf,name,'png')
close


figure
surf(X,Y,compare_lambda_2)
xlabel('gamma_1st=')
ylabel('gamma_2nd=')
zlabel('relative error of lambda_2')
name=[Bigname,'test_lambda_2'];
saveas(gcf,name,'png')
close
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


















