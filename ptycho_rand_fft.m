


function    Y =ptycho_rand_fft(X,os_rate,num_mask,mask,l_patch,x_c_p,y_c_p,bd_con)
   %%%% be aware x_c_p, y_c_p are matrix
  half_patch= (l_patch-1)/2;
  [Na,Nb]=size(X);
  
switch bd_con
    case 1 % periodic boundary condtion 
        
     Big_X=kron(ones(3,3),X);
     
     [subNa,subNb]=size(x_c_p);
     %%%% size of each patch in fourier space = os_rate* l_patch-by-
     %%%% os_rate*l_patch
     Y=zeros(os_rate*l_patch*subNa,os_rate*l_patch*subNb);
    for i = 1:subNb
        for j =1:subNa
            center_x = Na+ x_c_p(j,i);
            center_y = Nb+ y_c_p(j,i);
            piece=Big_X(center_x-half_patch:center_x+half_patch,center_y-half_patch:center_y+half_patch);
            piece_Y=Nos_fft(piece,os_rate,num_mask,mask);
            
            Y((j-1)*os_rate*l_patch+1 : j*os_rate*l_patch,  (i-1)*os_rate*l_patch+1 : i*os_rate*l_patch )=...
                piece_Y;
            
        end


    end
    
    
end










end