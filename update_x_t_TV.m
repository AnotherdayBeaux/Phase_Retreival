
%%%% the x_t inner loop in DR 1st step.




function  [x_t,count,test_cv]=update_x_t_TV(x_t, lambda, tau, tol, MaxIter)
    [Na, Nb]=size(x_t);
    p_n=zeros(Na,Nb,2);
    test_cv=1;
    Reb8p=x_t/lambda;
    count=1;
while test_cv > tol && count< MaxIter
    kk=grad_div(p_n,Reb8p);
    p_nplus1=(p_n+tau*kk)./(1+tau*abs(kk));
    test_cv=norm(div2(p_nplus1)-div2(p_n),'fro')/norm(div2(p_n),'fro');
    p_n=p_nplus1;
    count=count+1;
    
end

x_t=x_t-lambda*div2(p_n);

end