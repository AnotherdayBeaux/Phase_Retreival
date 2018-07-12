%%%%%% image print function for DR/ADMM method  
%%%%%%
%%%%%% '||A^Ax_t-y_t||_2/||y_t||_2'


% semi = 0; loglog plot    when qudratic g. 
% semi = 1; semilogy plot  when pois loglikelihood.
% showim = 1; show termial image
% showim = 0; not show terminal image


% Name = location



function image_print(Name,count,semi,residual_DR,dev,Lambda_2,x_t,coer)
switch semi
    case 1
    figure
    semilogy(1:count-1,residual_DR(1:count-1),'g')
    title('semily Relative error')
    xlabel('Iter')
    ylabel('|||A^{*}Ax|-b||_2/||b||_2')
    name=[Name,'/relative_error'];
    saveas(gcf,name,'png')
    close
    
    if nargin > 5
        figure
        semilogy(1:count-1,Lambda_2(1:count-1),'g')
        title('fixed point test')
        xlabel('Iter')
        ylabel('||AAx_t-x_t||/||x_t||')
        name=[Name,'/fixe_point_test'];
        saveas(gcf,name,'png')
        close
    end

    
    if nargin > 6 
        figure
        imshow(abs(x_t)*coer/255);
        title('recovered image')
        name=[Name,'/recovered_image'];
        saveas(gcf,name,'png')
        close
    end
    
    if nargin >4
        figure
        semilogy(1:count-1,dev(1:count-1),'g')
        title('semily deviation from Image')
        xlabel('Iter')
        ylabel('||IM-x_t||_2/||IM||_2')
        name=[Name,'/deviation from image'];
        saveas(gcf,name,'png')
        close
    end

    
    if count > 40000
        surname=[Name,'/after40000'];
        mkdir(surname)
        figure
        semilogy(count-20000:count-1,residual_DR(count-20000:count-1),'g')
        title('semily Relative error_last_20000')
        xlabel('Iter')
        ylabel('|||A^{*}Ax|-b||_2/||b||_2')
        name=[Name,'/relative_error_last_20000'];
        saveas(gcf,name,'png')
        close
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if nargin > 5
            figure
            semilogy(count-20000:count-1,Lambda_2(count-20000:count-1),'g')
            title('semily fixed_point_testlast20000')
            xlabel('Iter')
            ylabel('||AAx_t-x_t||/||x_t||')
            name=[Name,'/fixed_point_test_last20000'];
            saveas(gcf,name,'png')
            close
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure
        semilogy(count-20000:count-1,dev(count-20000:count-1),'g')
        title('semily deviation from Image')
        xlabel('Iter')
        ylabel('||IM-x_t||_2/||IM||_2')
        name=[Name,'/deviation_image_last20000'];
        saveas(gcf,name,'png')
        close

     
    end
    
    
    
    
    case 0       % loglog case
    figure
    loglog(1:count-1,residual_DR(1:count-1),'g')
    title(' Relative error')
    xlabel('Iter')
    ylabel('|||A^{*}Ax|-b||_2/||b||_2')
    name=[Name,'/relative_error'];
    saveas(gcf,name,'png')
    close
    
    if nargin > 5
        figure
        loglog(1:count-1,Lambda_2(1:count-1),'g')
        title('fixed point test')
        xlabel('Iter')
        ylabel('||AAx_t-x_t||/||x_t||')
        name=[Name,'/fixe_point_test'];
        saveas(gcf,name,'png')
        close
    end

    
    if nargin > 6 
        figure
        imshow(abs(x_t)*coer/255);
        title('recovered image')
        name=[Name,'/recovered_image'];
        saveas(gcf,name,'png')
        close
    end
        
    figure
    loglog(1:count-1,dev(1:count-1),'g')
    title(' deviation from Image')
    xlabel('Iter')
    ylabel('||IM-x_t||_2/||IM||_2')
    name=[Name,'/deviation from image'];
    saveas(gcf,name,'png')
    close

    
    if count > 40000
        surname=[Name,'/after40000'];
        mkdir(surname)
        figure
        loglog(count-20000:count-1,residual_DR(count-20000:count-1),'g')
        title(' Relative error last 20000')
        xlabel('Iter')
        ylabel('|||A^{*}Ax|-b||_2/||b||_2')
        name=[Name,'/relative_error_last_20000'];
        saveas(gcf,name,'png')
        close
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if nargin > 5
            figure
            loglog(count-20000:count-1,Lambda_2(count-20000:count-1),'g')
            title('fixed point testlast20000')
            xlabel('Iter')
            ylabel('||AAx_t-x_t||/||x_t||')
            name=[Name,'/fixed_point_test_last20000'];
            saveas(gcf,name,'png')
            close
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure
        loglog(count-20000:count-1,dev(count-20000:count-1),'g')
        title(' deviation from Image')
        xlabel('Iter')
        ylabel('||IM-x_t||_2/||IM||_2')
        name=[Name,'/deviation_image_last20000'];
        saveas(gcf,name,'png')
        close

     
    end
    



end