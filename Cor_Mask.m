% generate correlated mask convolve iid mask with square characteristic
% function 



function Cor_mask = Cor_Mask(m,c_l, phase_arg_input)


BigMask_phase=rand(m+c_l,m+c_l);
%zero_pad=c_l/2;
cl_2=floor(c_l/2);
BigMask_phase(cl_2+1: cl_2+m, cl_2+1: cl_2+m)=phase_arg_input;
%zero_ind=zeros(m+c_l,m+c_l);
%zero_ind(floor(zero_pad):m+c_l-floor(zero_pad)+1, floor(zero_pad):m+c_l-floor(zero_pad)+1)=1;

BigMask=exp((BigMask_phase-0.5)*2*pi*1i);

Cor_mask=zeros(m,m);
for i=1:c_l
    for j = 1:c_l
        Cor_mask=Cor_mask+BigMask(i:m+i-1,j:m+j-1);
         
        
    end
end

Cor_mask= Cor_mask./abs(Cor_mask);
Cor_mask=angle(Cor_mask)/2/pi;
end
