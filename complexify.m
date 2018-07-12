function compX=complexify(X)

[Na,X_Nb]=size(X);
Nb=X_Nb/2;
compX=X(:,1:Nb)+1i*X(:,Nb+1:end);


end