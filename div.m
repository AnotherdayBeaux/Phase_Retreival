%%%% calculate the divergence of p (chambolle total variation)


function div_p=div(p)


[Na,Nb,Nc]=size(p);
p_1=p(:,:,1);
p_2=p(:,:,2);
div_p=[p_1(1:Na-1,:);zeros(1,Nb)]-[zeros(1,Nb);p_1(2:Na,:)]+...
     [p_2(:,1:Nb-1),zeros(Na,1)]-[zeros(Na,1),p_2(:,2:Nb)];
                
                
                
                
end


function div_p=div_b(p)

[Na,Nb,Nc]=size(p);
p_1=p(:,:,1);
p_2=p(:,:,2);
div_p=[p_1(1:Na-1,:);zeros(1,Nb)]-[zeros(1,Nb);p_1(2:Na,:)]+...
     [p_2(:,1:Nb-1),zeros(Na,1)]-[zeros(Na,1),p_2(:,2:Nb)];
                
                
                
                
end
