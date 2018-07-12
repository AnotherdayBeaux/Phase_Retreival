%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Suppose A= a+b*i, deA returns deA=[a,b]. 
%%% 


function deA=decomplex(A)
     deA=[real(A),imag(A)];
end