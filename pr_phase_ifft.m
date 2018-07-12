
% pr_phase_fft applies inverse fourier transform to phase





function phase=pr_phase_ifft(fft_phase,os_rate,num_mask,X,l_patch,x_c_p,y_c_p)

[Na,Nb]=size(X);
Big_X= kron(ones(3,3), X);
half_patch = (l_patch-1)/2;
[subNa,subNb]=size(x_c_p);
phase=zeros(l_patch, l_patch);
for i =1: subNb
    for j=1: subNa
        center_x=x_c_p(j,i)+Na; center_y=y_c_p(j,i)+Nb;
        mask=Big_X(center_x-half_patch:center_x+half_patch, center_y-half_patch: center_y+half_patch);
        fft_patch_phase=fft_phase((j-1)*l_patch*os_rate+1:j*l_patch*os_rate, (i-1)*l_patch*os_rate+1:i*l_patch*os_rate);
        ifft_patch_phase=Nos_ifft(fft_patch_phase,os_rate,num_mask,mask);
        phase=phase+ifft_patch_phase;
        
    end
end
   


end