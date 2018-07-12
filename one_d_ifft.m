function IM= one_d_ifft(Y, os_rate, num_mask, mask)


    num_of_masks=floor(num_mask/2)+1;
    
    
switch num_mask
    
    case 0 
        
        [Na, Nb]= size(mask);
        X=zeros(Na*Nb,1);
        
        
        
        for i = 1: os_rate
            LL=diag(exp(2*(i-1)*1i*pi/os_rate * (0:1/(Na*Nb) : (Na*Nb-1)/Na/Nb)));
            X1=LL*ifft(Y(i:os_rate:end));
            X=X+X1;  
        end
        X=X*sqrt(Na*Nb);
        
        
    case 1
        [Na, Nb]= size(mask);
        X=zeros(Na*Nb,1);
        
        
        for i = 1: os_rate
            LL=diag(exp(2*(i-1)*1i*pi/os_rate * (0:1/(Na*Nb) : (Na*Nb-1)/Na/Nb)));
            X1=LL*ifft(Y(i:os_rate:end));
            X=X+X1;
        end
        X=conj(mask).*X;
        X=X*sqrt(Na*Nb);
           
    case 2
        [Na,Nb,Nc]=size(mask);
        X=zeros(Na*Nb,1);
        X2=zeros(Na*Nb,1);
        for i = 1: os_rate
            LL=diag(exp(2*(i-1)*1i*pi/os_rate * (0:1/(Na*Nb) : (Na*Nb-1)/Na/Nb)));
            X1=LL*ifft(Y(i:os_rate:Na*Nb*os_rate));
            X=X+X1;
            X22=LL*ifft(Y(Na*Nb*os_rate+i:os_rate:Na*Nb*os_rate*2));
            X2=X2+X22;
        end
        IM1=conj(mask(:,:,1)).*reshape(X,Na,Nb);
        IM2=conj(mask(:,:,2)).*reshape(X2,Na,Nb);
        IM=IM1+IM2;
        IM=IM*sqrt(Na*Nb); 
        
        
        
        
    case 3
        [Na,Nb,Nc]=size(mask);
        X=zeros(Na*Nb,1);
        X2=zeros(Na*Nb,1);
        for i = 1: os_rate
            LL=diag(exp(2*(i-1)*1i*pi/os_rate * (0:1/(Na*Nb) : (Na*Nb-1)/Na/Nb)));
            X1=LL*ifft(Y(i:os_rate:Na*Nb*os_rate));
            X=X+X1;
            X22=LL*ifft(Y(Na*Nb*os_rate+i:os_rate:Na*Nb*os_rate*2));
            X2=X2+X22;
        end
        IM1=conj(mask(:,:,1)).*reshape(X,Na,Nb);
        IM=IM1+reshape(X2,Na,Nb);
        IM=IM*sqrt(Na*Nb); 
              
end
        
IM=IM/sqrt(num_of_masks*os_rate);
end