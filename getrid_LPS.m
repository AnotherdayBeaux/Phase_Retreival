%%% this function minimize the following:

%   min_{k, theta} ||f_0 - exp(i theta) * exp(-i 2pi * k*n)f_i||_2

% the theta and k minimization can be done separately 





function k= getrid_LPS(f_0, f_k)

[Ia,Ib]=size(f_0);
A=(1:Ia)'*ones(1,Ib);
B=ones(1,Ia)'*(1:Ib);

exp_2pikn=f_0./f_k;           %%%%% exp(-2i pi k n)


%fftexp_2pikn=fft2(exp_2pikn);
%[k1,k2] = find(abs(fftexp_2pikn)== max(max(abs(fftexp_2pikn))));






exp_kn_matrix=exp_2pikn./exp(2i*pi*(k1*A+k2*B)/Ia);


kn_matrix=log(exp_kn_matrix)/(2i*pi) *Ia;


C=kn_matrix;


tr_AA= trace(A*A');  %a
tr_AB= trace(A*B');  %b
tr_AC= trace(A*C');  %c 

tr_BB= trace(B*B');  %f
tr_CB= trace(C*B');  %e
abbf=[tr_AA, tr_AB;tr_AB, tr_BB];
ce=[tr_AC; tr_CB];

k=inv(abbf)*ce;

k=k+[k1;k2];

%k=[k1,k2];

end
