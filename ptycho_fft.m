function    Y =ptycho_fft(X,os_rate,num_mask,mask,subim,overlap)
    [Na,Nb]=size(X);
     subinterv=floor(subim*overlap);
     subNa=Na/subinterv;
     subNb=Nb/subinterv;
     Y=zeros(os_rate*subim*(subNa-1),os_rate*subim*(subNb-1));
    for i = 1:subNb-1
        for j =1:subNa-1
            
            piece=X((j-1)*subinterv+1:(j+1)*subinterv,(i-1)*subinterv+1:(i+1)*subinterv);
            piece_Y=Nos_fft(piece,os_rate,num_mask,mask);
            
            Y((j-1)*os_rate*subim+1 : j*os_rate*subim,  (i-1)*os_rate*subim+1 : i*os_rate*subim )=...
                piece_Y;
            



        end


    end
end