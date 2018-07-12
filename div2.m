% calculate the discrete divergence of a function 3-d case




function div_p=div2(p)


[Na,Nb,Nc]=size(p);

div_p=zeros(Na,Nb);
div_p(2:Na-1,2:Nb-1)=p(2:Na-1,2:Nb-1,1)-p(1:Na-2,2:Nb-1,1)+p(2:Na-1,2:Nb-1,2)-p(2:Na-1,1:Nb-2,2);

div_p(1,1)=p(1,1,1)+p(1,1,2);
div_p(1,Nb)=p(1,Nb,1)-p(1,Nb-1,2);
div_p(Na,1)=-p(Na-1,1,1)+p(Na,1,2);
div_p(Na,Nb)=-p(Na-1,Nb,1)-p(Na,Nb-1,2);

div_p(1,2:Nb-1)=p(1,2:Nb-1,1)+p(1,2:Nb-1,2)-p(1,1:Nb-2,2);
div_p(Na,2:Nb-1)=-p(Na-1,2:Nb-1,1)+p(Na,2:Nb-1,2)-p(Na,1:Nb-2,2);
div_p(2:Na-1,1)=p(2:Na-1,1,1)-p(1:Na-2,1,1)+p(2:Na-1,1,2);
div_p(2:Na-1,Nb)=p(2:Na-1,Nb,1)-p(1:Na-2,Nb,1)-p(2:Na-1,Nb-1,2);



end