
phase_=ones(l_patch,l_patch);
fft_phase_=pr_phase_fft(phase_,os_rate,num_mask,IM,l_patch,x_c_p,y_c_p); 
nor_phase=pr_phase_ifft(fft_phase_,os_rate,num_mask,IM,l_patch,x_c_p,y_c_p);
    
est_phase=1/2*rand(l_patch,l_patch)-1/4;
mask= exp(2i*pi*(phase_arg));
mask_estimate= exp(2i*pi*(phase_arg+1*est_phase));
count_DR=1;  

nor_mask=mask_estimate.*sqrt(nor_phase);
fft_phase=pr_phase_fft(nor_mask./sqrt(nor_phase),os_rate,num_mask,IM,l_patch,x_c_p,y_c_p);
update_phase=1;    
    while  update_phase < 60
            nor_mask_tplus1=pr_phase_ifft(fft_phase,os_rate,num_mask,IM,l_patch,x_c_p,y_c_p)./sqrt(nor_phase);
            mask_tplus1=nor_mask_tplus1./abs(nor_mask_tplus1);
            nor_mask_tplus1= sqrt(nor_phase).*mask_tplus1;
            fft_phase_1=pr_phase_fft(nor_mask_tplus1./sqrt(nor_phase),os_rate,num_mask,IM,l_patch,x_c_p,y_c_p);
            
            Q_phase=2*fft_phase_1-fft_phase;
            Q=abs(Q_phase);
            rot_z=(Q/gamma+sqrt(Q.^2/gamma^2+8*(2+1/gamma)*(Z)))/(4+2/gamma);
            z_mask_tplus1 = rot_z.* Q_phase./Q;
            
            
            fft_phase=fft_phase+ z_mask_tplus1   -fft_phase_1;

            ee=norm(mask(:)'*mask_tplus1(:))/(mask(:)'*mask_tplus1(:));
            rel_mask=norm(ee*mask_tplus1 - mask,'fro')/norm(mask,'fro');
            fprintf('u_p=%d\n rel_mask=%f\n',update_phase,rel_mask);
            update_phase=update_phase+1;
            
    end