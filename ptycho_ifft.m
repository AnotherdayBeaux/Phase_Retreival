function    X =ptycho_ifft(Y,os_rate,num_mask,mask,subim,overlap)
            
  [Y_Na,Y_Nb]=size(Y);
  
  
  subinterv=floor(subim*overlap);
     
  Na=(Y_Na/os_rate/subim+1)*subinterv;
  Nb=(Y_Nb/os_rate/subim+1)*subinterv;
  
  subNa=Na/subinterv;
  subNb=Nb/subinterv;
  
  X=zeros(Na,Nb);
     for  i = 1:subNb-1
         for j= 1:subNa-1
             piece_Y=Y((j-1)*os_rate*subim+1 : j*os_rate*subim,  (i-1)*os_rate*subim+1 : i*os_rate*subim);
             piece_X=Nos_ifft(piece_Y,os_rate,num_mask,mask);
             X((j-1)*subinterv+1:(j+1)*subinterv,(i-1)*subinterv+1:(i+1)*subinterv)=...
                 X((j-1)*subinterv+1:(j+1)*subinterv,(i-1)*subinterv+1:(i+1)*subinterv)+piece_X;
             
             
         end
     end
     
     
     
end
  