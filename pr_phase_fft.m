
% pr_phase_fft applies fourier transform to phase
% pr_phase_fft(phase,os_rate,num_mask,X,l_patch,x_c_p,y_c_p)





function fft_phase=pr_phase_fft(phase,os_rate,num_mask,X,l_patch,x_c_p,y_c_p)


[Na,Nb]=size(X);
Big_X= kron(ones(3,3), X);
half_patch = (l_patch-1)/2;
[subNa,subNb]=size(x_c_p);
fft_phase=zeros(os_rate*l_patch*subNa,os_rate*l_patch*subNb);
for i =1: subNb
    for j=1: subNa
        center_x=x_c_p(j,i)+Na; center_y=y_c_p(j,i)+Nb;
        mask=Big_X(center_x-half_patch:center_x+half_patch, center_y-half_patch: center_y+half_patch);
        fft_patch_phase=Nos_fft(phase,os_rate,num_mask,mask);
        fft_phase(os_rate*l_patch*(j-1)+1: os_rate*l_patch*j, os_rate*l_patch*(i-1)+1: os_rate*l_patch*i )=fft_patch_phase;
        
    end
end
   


end