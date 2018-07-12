function obj=x_update_TV(xx,y_klambda,os_rate,num_mask,mask,tv,rho)  
               %%% xx= [real(x),imag(x)];

x=complexify(xx);
absx=abs(x);
[Na,Nb]=size(x);
c=zeros(Na,1);
d=zeros(1,Nb);
x_shift_horizon=[absx(:,2:end),c];
x_shift_vertical=[absx(2:end,:);d];
TV_1=x-x_shift_horizon;
TV_1=[TV_1(:,1:end-1),c];
TV_2=x-x_shift_vertical;
TV_2=[TV_2(1:end-1,:);d];

fx=os_fft(x,os_rate,num_mask,mask); %% complex form

% obj=norm(fx-y_klambda,'fro')^2+tv*sum(sum(abs(TV_1)+abs(TV_2)));    % l1 norm

obj=rho/2*norm(fx-y_klambda,'fro')^2+tv*sum(sum(sqrt(abs(TV_1).^2+abs(TV_2).^2)));    % l2 norm


end

%%%%%%%%%%%%%

