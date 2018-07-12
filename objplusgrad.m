function [obj,grad_obj] = objplusgrad(x,c1,c2a,c2b,zij)
   obj=c1*(x(1)^2+x(2)^2)+c2a*x(1)+c2b*x(2)-zij*log(x(1)^2+x(2)^2);
   grad_obj= [2*c1*x(1)+c2a-2*zij*x(1)/(x(1)^2+x(2)^2)
            2*c1*x(2)+c2b-2*zij*x(2)/(x(1)^2+x(2)^2)];
   
end




