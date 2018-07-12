%%% when the obj function in 2nd step is 4th order gaussian. 
%%% where p = (2Z+1)/2, Z= abs(Y)^2; q= -Q/2, Q = abs(ylambdak)

function rot_zz = fourth_gaussian(Z,Q,gamma)
   p=(1/gamma-2*Z)/2;
   q=-Q/(2*gamma);
   rot_z = zeros(length(p),3);
   rot_zz=zeros(length(p),1);

   w=(-1+sqrt(3)*1i)/2;
   rot_z(:,1)= (-q/2+sqrt((q/2).^2+(p/3).^3)).^(1/3)+(-q/2-sqrt((q/2).^2+(p/3).^3)).^(1/3);
   rot_z(:,2)= w*(-q/2+sqrt((q/2).^2+(p/3).^3)).^(1/3)+w^2*(-q/2-sqrt((q/2).^2+(p/3).^3)).^(1/3);
   rot_z(:,3)= w^2*(-q/2+sqrt((q/2).^2+(p/3).^3)).^(1/3)+w*(-q/2-sqrt((q/2).^2+(p/3).^3)).^(1/3);
   
   for i = 1: length(p)
       A=rot_z(i,:);
       
       obj = 1/2 * (abs(abs(A).^2-Z(i)).^2)+ 1/(2*gamma) * abs(abs(A)-Q).^2;
       a= find(obj==min(obj));   
       
       rot_zz(i)=rot_z(a(1));
          
   end
       
    
end