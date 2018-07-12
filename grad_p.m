
% periodic boundary 

function grad_u=grad_p(u)

[Na,Nb]=size(u);

grad_u=zeros(Na,Nb,2);
grad_u(:,:,1)=[u(2:Na,:);u(1,:)]-u;
grad_u(:,:,2)=[u(:,2:Nb),u(:,1)]-u;

end