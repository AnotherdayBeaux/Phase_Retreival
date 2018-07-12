%%%   solves the following optimization
%%%   f(z)+\rho/2 || z-lambda||^2
%%% 
%%%%  type_f = 1; f=|z|^2-blog(|z|^2)
%%%%  type_f = 0; f=1/2 |||z|-\sqrt{b}||^2


%%% poison likelihood.
%%% gaussi likelihood.

% J_F=(1+1/rho partial f)^-1



function z_kplus1=J_F(b,rho,lambda,RanS,type_f)
      Q=abs(lambda);
      
      
      switch type_f
          case 1
        % POISSON
        rot_z=rho/(4+2*rho)*(Q +sqrt(Q.^2+8*b*(2+rho)/rho^2));
          case 0 
        % GAUSSIAN
        rot_z= (sqrt(b)+rho*Q)/(1+rho);
      end
            
        z_kplus1=rot_z.*lambda./Q; 
        z_kplus1(isnan(z_kplus1))=0;
        z_kplus1=z_kplus1.*RanS;
        z_kplus1=z_kplus1.*RanS+abs(RanS-1).*lambda;

end