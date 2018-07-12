%%%% use FFT to calculate (A_star*D)*x 



function Y= Nos_fft(X,os_rate,num_mask,mask)
[Na,Nb]=size(X);
num_of_masks=floor(num_mask/2)+1;
switch num_mask   % 0  0 mask / Id mask
                  % 1  1 mask case
                  % 2  2 mask case 
                  % 3  1 and 1/2 mask case
                  % 4  3 mask, 1 id, 2 random ones.
    case 0
        Y=zeros(os_rate*Na,os_rate*Nb);
        %X=mask.*X;
        for i = 1: os_rate
            LL=diag(exp(-2*(i-1)*1i*pi/os_rate *(0:1/Na:(Na-1)/Na)));
            for j = 1:os_rate
                RR=diag(exp(-2*(j-1)*1i*pi/os_rate *(0:1/Nb:(Nb-1)/Nb)));
                X1=LL*X*RR;
                Y1=fft2(X1);
                Y(i:os_rate:os_rate*Na,j:os_rate:os_rate*Nb)=Y1;        
            end
        end
        
    case 1
        Y=zeros(os_rate*Na,os_rate*Nb);
        X=mask.*X;
        for i = 1: os_rate
            LL=diag(exp(-2*(i-1)*1i*pi/os_rate *(0:1/Na:(Na-1)/Na)));
            for j = 1:os_rate
                RR=diag(exp(-2*(j-1)*1i*pi/os_rate *(0:1/Nb:(Nb-1)/Nb)));
                X1=LL*X*RR;
                Y1=fft2(X1);
                Y(i:os_rate:os_rate*Na,j:os_rate:os_rate*Nb)=Y1;        
            end
        end
    
    case 2
        Y=zeros(os_rate*Na,os_rate*Nb*2);
        mask1X=mask(:,:,1).*X;
        mask2X=mask(:,:,2).*X;
        for i = 1: os_rate
            LL=diag(exp(-2*(i-1)*1i*pi/os_rate *(0:1/Na:(Na-1)/Na)));
            for j = 1:os_rate
                RR=diag(exp(-2*(j-1)*1i*pi/os_rate *(0:1/Nb:(Nb-1)/Nb)));
                X1=LL*mask1X*RR;
                X2=LL*mask2X*RR;
                Y1=fft2(X1);
                Y2=fft2(X2);
                Y(i:os_rate:os_rate*Na,j:os_rate:os_rate*Nb)=Y1;
                Y(i:os_rate:os_rate*Na,os_rate*Nb+j:os_rate:os_rate*Nb*2)=Y2;
            end
        end
        
        
        
    case 3
        Y=zeros(os_rate*Na,os_rate*Nb*2);
        mask1X=mask(:,:,1).*X;
        %mask2X=mask(:,:,2).*X;
        for i = 1: os_rate
            LL=diag(exp(-2*(i-1)*1i*pi/os_rate *(0:1/Na:(Na-1)/Na)));
            for j = 1:os_rate
                RR=diag(exp(-2*(j-1)*1i*pi/os_rate *(0:1/Nb:(Nb-1)/Nb)));
                X1=LL*mask1X*RR;
                X2=LL*X*RR;
                Y1=fft2(X1);
                Y2=fft2(X2);
                Y(i:os_rate:os_rate*Na,j:os_rate:os_rate*Nb)=Y1;
                Y(i:os_rate:os_rate*Na,os_rate*Nb+j:os_rate:os_rate*Nb*2)=Y2;
            end
        end
        
        
        
        case 4     %%%%%% 3 masks
        Y=zeros(os_rate*Na,os_rate*Nb*2);
        mask1X=mask(:,:,1).*X;
        mask2X=mask(:,:,2).*X;
        for i = 1: os_rate
            LL=diag(exp(-2*(i-1)*1i*pi/os_rate *(0:1/Na:(Na-1)/Na)));
            for j = 1:os_rate
                RR=diag(exp(-2*(j-1)*1i*pi/os_rate *(0:1/Nb:(Nb-1)/Nb)));
                X1=LL*mask1X*RR;
                X2=LL*mask2X*RR;
                X3=LL*X*RR;
                Y1=fft2(X1);
                Y2=fft2(X2);
                Y3=fft2(X3);
                Y(i:os_rate:os_rate*Na,j:os_rate:os_rate*Nb)=Y1;
                Y(i:os_rate:os_rate*Na,os_rate*Nb+j:os_rate:os_rate*Nb*2)=Y2;
                Y(i:os_rate:os_rate*Na,2*os_rate*Nb+j:os_rate:os_rate*Nb*3)=Y3;
            end
        end
        
     
      
end
        
        
        
        
Y=Y*1/sqrt(Na*Nb)/sqrt(num_of_masks*os_rate^2);

end
