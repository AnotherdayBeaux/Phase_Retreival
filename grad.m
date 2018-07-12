function grad_u=grad(u)

[Na,Nb]=size(u);

grad_u=zeros(Na,Nb,2);
grad_u(1:Na-1,:,1)=u(2:Na,:)-u(1:Na-1,:);
grad_u(:,1:Nb-1,2)=u(:,2:Nb)-u(:,1:Nb-1);

end