

% p_c = probe caliberation 
%



function    X =ptycho_rand_ifft(Y,Na,Nb,os_rate,num_mask,mask,l_patch,x_c_p,y_c_p, p_c)
  half_patch= (l_patch-1)/2;
  Big_X=zeros(3*Na, 3*Nb);
 
  
  
     
  
  [subNa, subNb]= size(x_c_p);
  J_vec=[subNb:-1:1;1:subNb];
  switch p_c
      case 0  % caliberation off
     for  i = 1:subNa
         rev= mod(i,2);
         for j = J_vec(rev+1,:)
             center_x=x_c_p(i,j)+Na; center_y=y_c_p(i,j)+Nb;
             piece_Y=Y((i-1)*os_rate*l_patch+1 : i*os_rate*l_patch,  (j-1)*os_rate*l_patch+1 : j*os_rate*l_patch);
             piece_X=Nos_ifft(piece_Y,os_rate,num_mask,mask);
             Big_X(center_x-half_patch:center_x+half_patch, center_y-half_patch: center_y+half_patch)=...
                 Big_X(center_x-half_patch:center_x+half_patch, center_y-half_patch: center_y+half_patch)+piece_X;
             
             
         end
     end
     
      case 1 % on
    
     piece_Y=Y(1 : os_rate*l_patch,  1 : os_rate*l_patch);
     piece_X_pre=Nos_ifft(piece_Y,os_rate,num_mask,mask);
     center_x_pre=x_c_p(1,1)+Na;
     center_y_pre=y_c_p(1,1)+Nb;
     for  i = 1:subNa
         
         rev= mod(i,2);
         for j = J_vec(rev+1,:)
             
             piece_Y=Y((i-1)*os_rate*l_patch+1 : i*os_rate*l_patch,  (j-1)*os_rate*l_patch+1 : j*os_rate*l_patch);
             piece_X=Nos_ifft(piece_Y,os_rate,num_mask,mask);
             center_x_start=x_c_p(i,j)+Na; center_y_start=y_c_p(i,j)+Nb;
             COS=zeros(11,11);
             s_x=0;
             for search_x= -5:5
                 s_x=s_x+1;
                 s_y=0;
                 for search_y = -5:5
                     s_y=s_y+1;
                    Big_X_olpre=zeros(3*Na, 3*Nb);
                    Big_X_olpost=zeros(3*Na, 3*Nb);
                    Big_X_olpre(center_x_pre-half_patch:center_x_pre+half_patch,center_y_pre-half_patch:center_y_pre+half_patch)=...
                    ones(l_patch,l_patch);
                    
                    
                    Big_X_olpost(center_x_start+search_x-half_patch:center_x_start+search_x+half_patch,center_y_start+search_y-half_patch:center_y_start+search_y+half_patch)=...
                    ones(l_patch,l_patch);
                    Big_X_ind=Big_X_olpre.*Big_X_olpost;

                    Big_X_pre=zeros(3*Na, 3*Nb);
                    Big_X_post=zeros(3*Na, 3*Nb);
                    Big_X_pre(center_x_pre-half_patch:center_x_pre+half_patch,center_y_pre-half_patch:center_y_pre+half_patch)=...
                    piece_X_pre;
                    Big_X_post(center_x_start+search_x-half_patch:center_x_start+search_x+half_patch,center_y_start+search_y-half_patch:center_y_start+search_y+half_patch)=...
                    piece_X;
                    cos=norm(Big_X_pre.*conj(Big_X_post),'fro')^2/norm(Big_X_ind.*Big_X_pre,'fro')/norm(Big_X_ind.*Big_X_post,'fro');
                    COS(s_x,s_y)=cos;  
                 end
             end
             k=max(COS(:));
             [a,b]=find(COS==k);
             center_x=center_x_start-6+a(1);
             center_y=center_y_start-6+b(1);

             
             
             Big_X(center_x-half_patch:center_x+half_patch, center_y-half_patch: center_y+half_patch)=...
                 Big_X(center_x-half_patch:center_x+half_patch, center_y-half_patch: center_y+half_patch)+piece_X;
             
             piece_X_pre=piece_X;
             center_x_pre=center_x;
             center_y_pre=center_y;
         end
     end  
     
     
     
     
     
  end

  X=zeros(Na,Nb);
  for i = 1:3
      for j=1:3
        X=Big_X((i-1)*Na+1:i*Na,(j-1)*Nb+1:j*Nb)+X;
      end
  end
  
  
     
     
     
end