%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                  %%%
%%%       sector constrains          %%%
%%% x_k proj  [-alpha \pi, beta \pi] %%%
%%%                                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                    %
%%%  case 0 general sector constrains; %
%%%       1 real  positivity;          %
%%%       2 complex positivity;        % 
%%%       3 non constrains;            %
%%%                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function projed_x_k=projection(x_k,sector_alpha,sector_beta,sector_mode)
if nargin <2
    sector_mode = 0;
end
projed_x_k=zeros(size(x_k));

switch sector_mode
    case 0       
        [Na,Nb]=size(x_k);
        for i = 1: Na
            for j= 1:Nb
                phase_x_k_ij=phase(x_k(i,j));
                if     sector_beta*pi >= phase_x_k_ij  && phase_x_k_ij >= -sector_alpha*pi
                    projed_x_k(i,j)=x_k(i,j);
                elseif  (sector_beta+0.5)*pi >= phase_x_k_ij   && phase_x_k_ij >= sector_beta*pi
                    projed_x_k(i,j)=real(x_k(i,j)*exp(-sector_beta*pi*1i))*exp(sector_beta*pi*1i);
                elseif -sector_alpha*pi >= phase_x_k_ij  && phase_x_k_ij >= -(sector_alpha+0.5)*pi
                    projed_x_k(i,j)=real(x_k(i,j)*exp(sector_alpha*pi*1i))*exp(-sector_alpha*pi*1i);
                else
                    projed_x_k(i,j)=0;
                end           
            end  
        end
        
        
    case 1
        projed_x_k=max(real(x_k),0); 
        
    case 2
        projed_x_k=max(real(x_k),0)+1i*max(imag(x_k),0);
        
    case 3
        projed_x_k=x_k;
        
        
end
        


end