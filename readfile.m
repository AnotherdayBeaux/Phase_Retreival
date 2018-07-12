
bigname='/Users/Beaux/Downloads/fig_nonoise_Cameraman_compare_ga1_ga2';
kk=dir('/Users/Beaux/Downloads/fig_nonoise_Cameraman_compare_ga1_ga2');


iter=zeros(25,149);
for i = 4:3728
    foldername=kk(i).name;
    textname=[bigname,'/',foldername,'/configration.txt'];
    [a1,a2]=textread(textname,'%s%f','delimiter','=');
    
    gamma_1=a2(1);
    
    gamma_2=a2(2);
    
    iter(round(gamma_1-5),round((gamma_2-0.2)/0.2))=a2(4);
    
    
    
    
end

figure
ga1=6:30;
ga1=1./ga1;
ga2=0.4:0.2:30;
[X,Y]=meshgrid(ga1,ga2);
surf(ga1,ga2,iter');
xlabel('gamma1 1./(6:30)')
ylabel('gamma2 0.4:0.2:30')
title('compare gamma1 gamma2 in pois')


figure
ga2=0.8:0.2:30;
[X,Y]=meshgrid(ga1,ga2);
surf(ga1,ga2,iter(:,3:end)');
xlabel('gamma1 1./(6:30)')
ylabel('gamma2 0.8:0.2:30')
title('compare gamma1 gamma2 in pois')


    