function div_p=div_p(p)

[Na,Nb,Nc]=size(p);
p_1=p(:,:,1);
p_2=p(:,:,2);
div_p=p_1-[p_1(Na,:);p_1(1:Na-1,:)]+...
     p_2-[p_2(:,Nb),p_2(:,1:Nb-1)];
                              
                
                
end
