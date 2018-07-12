%%% solve alpha *|| p ||_1 + rho/2 * || p - lambda||^2





function p_kplus1 = J_G(rho,alpha,lambda,type_l1)
            thresh=alpha/rho; 
            
            
            switch type_l1
                case 1
                p_kplus1=max(0,abs(lambda)-thresh).*lambda./abs(lambda);
                case 0
                re_inter=real(lambda);
                im_inter=imag(lambda);
            
                re_inter(re_inter>0)=1;
                re_inter(re_inter<0)=-1;
            
                im_inter(im_inter>0)=1;
                im_inter(im_inter<0)=-1;  
            
                real_inter=real(lambda);
                imag_inter=imag(lambda);
            
                p_kplus1=max(0,abs(real_inter)-thresh).*re_inter+1i*max(0,abs(imag_inter)-thresh).*im_inter;
                
                
            end
end


