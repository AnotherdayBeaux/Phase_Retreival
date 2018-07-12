%%%% this method get the convergence rate of last part of iteration.


function cr=get_cr(seq,last)


subseq=seq(end-2*last+1:end);
cr=(sum(subseq(1:last)-subseq(last+1:end)))/last/last;
cr=exp(-abs(cr));
end