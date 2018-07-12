%%% x, y_k, c are Na-by-Nb-by-2 matrix,
%%% objective function of ||y-Ax||^2  
%%% want to use quadrag or fminunc to solve


function [obj,grad]=x_update(xx,y_klambda,os_rate,num_mask,mask)  
               %%% xx= [real(x),imag(x)];

x=complexify(xx);

fx=os_fft(x,os_rate,num_mask,mask); %% complex form

obj=norm(fx-y_klambda,'fro')^2;
        
xxx=2*os_rate^2 *x*(floor(num_mask/2)+1);
yyy=2*os_ifft(y_klambda,os_rate,num_mask,mask);
grad=decomplex(xxx-yyy);
grad=grad(:);



end


