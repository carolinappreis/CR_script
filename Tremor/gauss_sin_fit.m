clear all
cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data')
load('smooth_arc3.mat')

cr=smo_s(1,:)';
bar(cr)
c=cr;


x = cr;
y = normpdf(x,0,1);

plot(x,y)



clearvars -except cr 
pd_kernel = fitdist(cr,'Kernel')
pd_normal = fitdist(cr,'Normal')
% x=1:12;

 x = min(cr):0.1:max(cr);
pdf_kernel = pdf(pd_kernel,x);
pdf_normal = pdf(pd_normal,x);

plot(x,pdf_kernel,'Color','b','LineWidth',2);
hold on;
plot(x,pdf_normal,'Color','r','LineStyle',':','LineWidth',2);
legend('Kernel Distribution','Normal Distribution','Location','SouthEast');
