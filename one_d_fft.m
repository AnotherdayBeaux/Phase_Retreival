function Y= one_d_fft(IM, os_rate, num_mask, mask)
    [Na,Nb]=size(IM);
    IM_vec=IM(:);
    num_of_masks=floor(num_mask/2)+1;
    
    
switch num_mask
    
    case 0 
        
        Y=zeros(os_rate*Na*Nb,1);
        
        
        for i = 1: os_rate
            LL=diag(exp(-2*(i-1)*1i*pi/os_rate * (0:1/(Na*Nb) : (Na*Nb-1)/Na/Nb)));
            X1=LL*IM_vec;
            Y1=fft(X1);
            Y(i:os_rate:(Na*Nb*os_rate))=Y1;   
        end
        
        
        
    case 1
        
        Y=zeros(os_rate*Na*Nb,1);
        
        IM_vec=mask(:).*IM_vec;
        
        for i = 1: os_rate
            LL=diag(exp(-2*(i-1)*1i*pi/os_rate * (0:1/(Na*Nb) : (Na*Nb-1)/Na/Nb)));
            X1=LL*IM_vec;
            Y1=fft(X1);
            Y(i:os_rate:(Na*Nb*os_rate))=Y1;   
        end
        
        
        
        
    case 2
        Y=zeros(os_rate*Na*Nb *2,1);
        mask1X=mask(:,:,1).*IM;
        mask2X=mask(:,:,2).*IM;
        IM_vec1=mask1X(:);
        IM_vec2=mask2X(:);
        
        for i = 1: os_rate
            LL=diag(exp(-2*(i-1)*1i*pi/os_rate * (0:1/(Na*Nb) : (Na*Nb-1)/Na/Nb)));
            X1=LL*IM_vec1;
            X2=LL*IM_vec2;
            Y1=fft(X1);
            Y2=fft(X2);
            Y(i:os_rate:(Na*Nb*os_rate))=Y1; 
            Y(Na*Nb*os_rate+i : os_rate:Na*Nb*os_rate*2 )= Y2;
        end
        
        
        
        
    case 3
        Y=zeros(os_rate*Na*Nb *2,1);
        mask1X=mask(:,:,1).*IM;
        IM_vec1=mask1X(:);
        
        for i = 1: os_rate
            LL=diag(exp(-2*(i-1)*1i*pi/os_rate * (0:1/(Na*Nb) : (Na*Nb-1)/Na/Nb)));
            X1=LL*IM_vec1;
            X2=LL*IM_vec;
            Y1=fft(X1);
            Y2=fft(X2);
            Y(i:os_rate:(Na*Nb*os_rate))=Y1; 
            Y(Na*Nb*os_rate+i : os_rate:Na*Nb*os_rate*2 )= Y2;
        end
        
        
        
        
        
        
        
end
        
        

Y=Y*1/sqrt(Na*Nb)/sqrt(num_of_masks*os_rate);
end







