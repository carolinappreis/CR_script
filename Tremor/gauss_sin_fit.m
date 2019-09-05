clear all
% cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data')
 cd('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data')
load('smooth_arc3.mat')
for i=1:10;
xdata = (1:12);
y=smo_s(i,:);
% bar(cr)
fs=fittype('a*sin(b*x)');
gs=fittype('a1*exp(-((x-b1)/c1).^2)');
s = fit(xdata',y',fs,'StartPoint',[max(y) 0.5]);
g = fit(xdata',y',gs,'StartPoint',[max(y) mean(xdata) 0.5*(max(xdata)-min(xdata))]);
figure(1)
subplot(10,1,i)
plot(s,xdata,y)
hold on 
yline(0)
legend('off')
ylim([-0.5 0.5])
box('off')

figure(2)
subplot(10,1,i)
plot(g,xdata,y)
hold on 
yline(0)
legend('off')
ylim([-0.5 0.5])
box('off')
clear y g s 
end

% cl

bar plot(smo_s(1,:))
data=smo_s(1,:);
zd=zscore(data);

xd=1:12;
s=max(data)*sin(0.5*xd+2*pi);
zs=zscore(s);

g = gaussmf(xd,[3 6]);
zg=zscore(g);

plot(xd,zd)
hold on
plot(xd,zg)
plot(xd,zs)

corr(zd,zg)







% % % for i=1:10;
% % % subplot(10,1,i)
% % % data=smo_s(i,:);
% % % xdata=1:12;
% % % fc = sin(xdata);
% % % %
% % % ft = fittype('a*sin(b*x)');
% % % f = fit(data',fc',ft,'StartPoint',[1 1]);
% % % hold on
% % % plot(f,data,fc)
% % % legend('off')
% % % end





tol=1e-6;

M=length(pha);

%gaussian param

A=max(data);
x0=mean(xdata);
s=0.5*(max(xdata)-min(xdata));
f=A*exp(-((x0)/s).^2);
plot(xdata,f)










err=inf;

while err>tol
    
    f=A*exp(-((cr-x0)/s).^2);
    
    d=cr-f;
    
    %construct matriz Z
   
    z1=f/A;
    z2=2*f.*(pha-x0)/s^2;
    z3=2*f.*(pha-x0).^2/s^3;
    z=[z1 z2 z3];
    
    %update coeficient 
    da=(z.'*z)\(z.'*d);
    A=A+da(1);
    x0=x0+da(2);
    s=s+da(3);
    
    %calculate error
    
    err=sum(abs(da./[A;x0;s]));
end

figure()
x=linspace(0,4,1000);
f=A*exp(-((x-x0)/s).^2);
plot(x,f)
hold
plot(pha,cr,'xr')



%
%
% x = cr;
% y = normpdf(x,0,1);
%
% plot(x,y)
%
%
%
% clearvars -except cr
% pd_kernel = fitdist(cr,'Kernel')
% pd_normal = fitdist(cr,'Normal')
% % x=1:12;
%
%  x = min(cr):0.1:max(cr);
% pdf_kernel = pdf(pd_kernel,x);
% pdf_normal = pdf(pd_normal,x);
%
% plot(x,pdf_kernel,'Color','b','LineWidth',2);
% hold on;
% plot(x,pdf_normal,'Color','r','LineStyle',':','LineWidth',2);
% legend('Kernel Distribution','Normal Distribution','Location','SouthEast');
