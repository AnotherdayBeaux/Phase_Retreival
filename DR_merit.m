%%%% test the convergence of merit function

function D=DR_merit(y_t,z_t,x_t)

D = 1/2*norm((abs(y_t)-b),'fro')^2 + 0+ ...
    1/2/gamma *(norm(x_t-y_t,'fro')^2-norm(x_t-z_t,'fro')^2);

end









