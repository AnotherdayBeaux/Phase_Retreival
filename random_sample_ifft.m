%%% random sample fourier inverse transform. IFFT(RanS*y), RanS is a diag
%%% matrix with 0,1 on diag. 



function X=random_sample_ifft(Y,os_rate,num_mask,mask,RanS)
% 0  0 mask / Id mask
% 1  1 mask case
% 2  2 mask case 
% 3  1 and 1/2 mask case
% 4  1 id, 2 random ones.
[Y_Na,Y_Nb]=size(Y);
Na=Y_Na/os_rate;
Y=RanS.*Y;
num_of_masks=floor(num_mask/2)+1;

switch num_mask
    case 0      % 0  0 mask / Id mask  
Nb=Y_Nb/os_rate;
X=zeros(Na,Nb);
for i = 1:os_rate
    LL=diag(exp(2*(i-1)*1i*pi/os_rate *(0:1/Na:(Na-1)/Na)));
    for j = 1:os_rate
        RR=diag(exp(2*(j-1)*1i*pi/os_rate *(0:1/Nb:(Nb-1)/Nb)));
        X1=LL*ifft2(Y(i:os_rate:Y_Na,j:os_rate:Y_Nb))*RR;
        X=X+X1;
    end
end
X=X*sqrt(Na*Nb);


    case 1   % 1  1 mask case
        Nb=Y_Nb/os_rate;
        X=zeros(Na,Nb);
for i = 1:os_rate
    LL=diag(exp(2*(i-1)*1i*pi/os_rate *(0:1/Na:(Na-1)/Na)));
    for j = 1:os_rate
        RR=diag(exp(2*(j-1)*1i*pi/os_rate *(0:1/Nb:(Nb-1)/Nb)));
        X1=LL*ifft2(Y(i:os_rate:Y_Na,j:os_rate:Y_Nb))*RR;
        X=X+X1;
    end
end
X=conj(mask).*X;
X=X*sqrt(Na*Nb);


    case 2   % 2  2 mask case 
        Nb=Y_Nb/os_rate/2;
        X=zeros(Na,Nb);
        X2=zeros(Na,Nb);
for i = 1:os_rate
    LL=diag(exp(2*(i-1)*1i*pi/os_rate *(0:1/Na:(Na-1)/Na)));
    for j = 1:os_rate
        RR=diag(exp(2*(j-1)*1i*pi/os_rate *(0:1/Nb:(Nb-1)/Nb)));
        X1=LL*ifft2(Y(i:os_rate:Y_Na,j:os_rate:Y_Nb/2))*RR;
        X=X+X1;
        X22=LL*ifft2(Y(i:os_rate:Y_Na,Y_Nb/2+j:os_rate:Y_Nb))*RR;
        X2=X2+X22;
    end
end
X=conj(mask(:,:,1)).*X;
X2=conj(mask(:,:,2)).*X2;
X=X+X2;
X=X*sqrt(Na*Nb); 



    case 3  % 1 1/2 mask  mask(:,:,1)=rand  mask(:,:,2)=Id
        Nb=Y_Nb/os_rate/2;
        X=zeros(Na,Nb);
        X2=zeros(Na,Nb);
for i = 1:os_rate
    LL=diag(exp(2*(i-1)*1i*pi/os_rate *(0:1/Na:(Na-1)/Na)));
    for j = 1:os_rate
        RR=diag(exp(2*(j-1)*1i*pi/os_rate *(0:1/Nb:(Nb-1)/Nb)));
        X1=LL*ifft2(Y(i:os_rate:Y_Na,j:os_rate:Y_Nb/2))*RR;
        X=X+X1;
        X22=LL*ifft2(Y(i:os_rate:Y_Na,Y_Nb/2+j:os_rate:Y_Nb))*RR;
        X2=X2+X22;
    end
end
X=conj(mask(:,:,1)).*X;
X=X+X2;
X=X*sqrt(Na*Nb); 


    case 4   % 4  1 id, 2 random ones.
        Nb=Y_Nb/os_rate/num_of_masks;
        X=zeros(Na,Nb);
        X2=zeros(Na,Nb);
        X3=zeros(Na,Nb);
for i = 1:os_rate
    LL=diag(exp(2*(i-1)*1i*pi/os_rate *(0:1/Na:(Na-1)/Na)));
    for j = 1:os_rate
        RR=diag(exp(2*(j-1)*1i*pi/os_rate *(0:1/Nb:(Nb-1)/Nb)));
        X1=LL*ifft2(Y(i:os_rate:Y_Na,j:os_rate:Y_Nb/3))*RR;
        X=X+X1;
        X22=LL*ifft2(Y(i:os_rate:Y_Na,Y_Nb/3+j:os_rate:(Y_Nb*2/3)))*RR;
        X2=X2+X22;
        X33=LL*ifft2(Y(i:os_rate:Y_Na,(Y_Nb*2/3)+j:os_rate:Y_Nb))*RR;
        X3=X3+X33;
    end
end
X=conj(mask(:,:,1)).*X;
X2=conj(mask(:,:,2)).*X2;
X=X+X2+X3;
X=X*sqrt(Na*Nb); 
               
end
        
X=X/sqrt(num_of_masks*os_rate^2);

end