% project x_t to periodic space



function p_IM= proj_periodic(IM,subim)
   
   [Na,Nb]=size(IM);
   % A B C D
  
   p_IM(subim+1:Na-subim,1:subim)=(IM(subim+1:Na-subim,1:subim)+IM(subim+1:Na-subim,Nb-subim+1:Nb))/2;
   p_IM(subim+1:Na-subim,Nb-subim+1:Nb)=(IM(subim+1:Na-subim,1:subim)+IM(subim+1:Na-subim,Nb-subim+1:Nb))/2;
   
   p_IM(1:subim,subim+1:Nb-subim)=(IM(1:subim,subim+1:Nb-subim)+IM(Na-subim+1:Na,subim+1:Nb-subim))/2;
   p_IM(Na-subim+1:Na,subim+1:Nb-subim)=(IM(1:subim,subim+1:Nb-subim)+IM(Na-subim+1:Na,subim+1:Nb-subim))/2;
   
   
   % E F H I 
   
   EFHI= (IM(1:subim,1:subim)+IM(1:subim,Nb-subim+1:Nb)+IM(Na-subim+1:Na,1:subim)+IM(Na-subim+1:Na,Nb-subim+1:Nb))/4;
   p_IM(1:subim,1:subim)=EFHI;
   p_IM(1:subim,Nb-subim+1:Nb)=EFHI;
   p_IM(Na-subim+1:Na,1:subim)=EFHI;
   p_IM(Na-subim+1:Na,Nb-subim+1:Nb)=EFHI;
   
   
   
   % rest
   p_IM(subim+1:Na-subim,subim+1:Nb-subim)=IM(subim+1:Na-subim,subim+1:Nb-subim);
   

end