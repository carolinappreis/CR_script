clear all
close all

load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\cleaned_rc12.mat')
load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\newnonstim2.mat')
load('C:\Users\creis\Documents\GitHub\CR_script\colour_pal.mat','blushred','squash');
%
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/cleaned_rc12.mat')
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/newnonstim2.mat')
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/labels_shift.mat')
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/arc_mediansplit_0120.mat')
% load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','squash');
cl=blushred;

%attention to number 13 and 14 of tt1 (they are the same pateint, 13 5
%trials with 5 pulses stim and 14, 5 trials with 1 pulse stim


%smoothed ARC'S
% iiii=[1 2 3 4 5 8 10 11 12 13 16 17 19 20];

cohort=[2 3 4 5 8 10 11 13 16 17];
dum=intersect(iiii,cohort);

pt=[];
for i=1:length(dum)
    pt=[pt find(iiii==dum(i))];
end

main=[1 1 3 1 3 3 3 3 1 1];
%   main=[1 1 1 1 1 1 1 1 1 1];
for pp=1:size(pt,2)
    for kk=1:size(tt1,2)
        
        curve=nanmedian(tt1{pt(pp),kk},1);
        SO=repmat(curve,1,3);
        SOUT=repmat((squeeze(nostimout(pt(pp),kk,:)))',1,3);
        for i=size(tt1{pt(pp),kk},2)+1:size(tt1{pt(pp),kk},2)*2
            smooth_c(1,i-12)=sum(SO(1,(i-1:i+1)))./length(SO(1,(i-1:i+1)));
            smooth_ns(1,i-12)=sum(SOUT(1,(i-1:i+1)))./length(SOUT(1,(i-1:i+1)));
        end
        
        
        smoo_all(pp,kk,:)=smooth_c;  clear smooth_c
        nsmoo_all(pp,kk,:)=smooth_ns; clear smooth_ns
        
    end
    
    
    
    smoo_main(pp,:)=squeeze(smoo_all(pp,1,:));
    nsmoo_main(pp,:)=squeeze(nsmoo_all(pp,1,:));
    s_main(pp,:)=nanmedian(tt1{pt(pp),1},1);
    ns_main(pp,:)=squeeze(nostimout(pt(pp),1,:));
end


%

for i=1:size(s_main,1);
    y= s_main(i,:);
    clearvars -except i y rs_gauss rs_sin cl s_main rs_lin
    f1=figure(1)
    subplot(1,10,i)
    bar(y,'FaceColor',cl,'EdgeColor',cl)
    hold on
    rsg=gauss_fit(y);
    rs_gauss(i,:)=[rsg.sse rsg.dfe];
    %     legend( 'ARC', 'gaussian fit', 'Location', 'NorthEast', 'Interpreter', 'none');
    %     legend('boxoff')
    %    xlabel( 'Stim phase', 'Interpreter', 'none' );
    %    ylabel( 'Amplitude change', 'Interpreter', 'none' );
    title(['r^2: ',num2str(rsg.rsquare)])
    %   ylim([-(max(abs(y)))-0.05 (max(abs(y))+0.05)])
    xlabel('');ylabel('');legend('off')
    box('off')
    
    f2=figure(2)
    subplot(1,10,i)
    bar(y,'FaceColor',cl,'EdgeColor',cl)
    hold on
    rss=sin_fit(y);
    rs_sin(i,:)=[rss.sse rss.dfe];
    xlabel('');ylabel('');legend('off')
    %     legend( 'ARC','Sin fit', 'Location', 'NorthEast', 'Interpreter', 'none');
    %     legend('boxoff')
    %     xlabel( 'Stim phase', 'Interpreter', 'none' );
    %     ylabel( 'Amplitude change', 'Interpreter', 'none' );
    title(['r^2: ',num2str(rss.rsquare)])
    %     ylim([-(max(abs(y)))-0.05 (max(abs(y))+0.05)])
    box('off')
    hold off
    
    
    f3=figure(3)
    subplot(1,10,i)
    bar(y,'FaceColor',cl,'EdgeColor',cl)
    hold on

    %    legend( 'ARC','Sin fit', 'Location', 'NorthEast', 'Interpreter', 'none');
    %    legend('boxoff')
    %    xlabel( 'Stim phase', 'Interpreter', 'none' );
    %    ylabel( 'Amplitude change', 'Interpreter', 'none' );
    rsl=poli_fit(y);
    rs_lin(i,:)=[rsl.sse rsl.dfe];
    xlabel('');ylabel('');legend('off')
    %     legend( 'ARC','Sin fit', 'Location', 'NorthEast', 'Interpreter', 'none');
    %     legend('boxoff')
    %     xlabel( 'Stim phase', 'Interpreter', 'none' );
    %     ylabel( 'Amplitude change', 'Interpreter', 'none' );
      title(['r^2: ',num2str(rsl.rsquare)])
    %     ylim([-(max(abs(y)))-0.05 (max(abs(y))+0.05)])
    box('off')
    hold off
    
end

for i=1:10
    n=12
    r1=rs_sin(i,1);
    r2=rs_lin(i,1);
    df1=rs_gauss(i,2);
    df2=rs_lin(i,2);
%     ctl_val=2.854;
    F=((r1-r2)./(df2-df1))./((r2)./(n-df2));
    fpdf(F,11,10)   
end
% for i=1:10
% y=smo_s(i,:);
% % bar(cr)
% fs=fittype('a*sin(b*x)+c');
% gs=fittype('a1*exp(-((x-b1)/c1).^2)');
% s = fit(xdata',y',fs,'StartPoint',[max(y) 0.5 0]);
% g = fit(xdata',y',gs,'StartPoint',[max(y) mean(xdata) 0.5*(max(xdata)-min(xdata))]);
% figure(1)
% subplot(10,1,i)
% plot(s,xdata,y)
% hold on
% yline(0)
% legend('off')
% ylim([-0.5 0.5])
% box('off')
%
% figure(2)
% subplot(10,1,i)
% plot(g,xdata,y)
% hold on
% yline(0)
% legend('off')
% ylim([-0.5 0.5])
% box('off')
% clear y g s
% end
%
% % cl
%
% bar plot(smo_s(1,:))
% data=smo_s(1,:);
% zd=zscore(data);
%
% xd=1:12;
% s=max(data)*sin(0.5*xd+2*pi);
% zs=zscore(s);
%
% g = gaussmf(xd,[3 6]);
% zg=zscore(g);
%
% plot(xd,zd)
% hold on
% plot(xd,zg)
% plot(xd,zs)
%
% corr(zd,zg)
%
%
%
%
%
%
%
% % % % for i=1:10;
% % % % subplot(10,1,i)
% % % % data=smo_s(i,:);
% % % % xdata=1:12;
% % % % fc = sin(xdata);
% % % % %
% % % % ft = fittype('a*sin(b*x)');
% % % % f = fit(data',fc',ft,'StartPoint',[1 1]);
% % % % hold on
% % % % plot(f,data,fc)
% % % % legend('off')
% % % % end
%
%
%
%
% tol=1e-6;
%
% M=length(pha);
%
% %gaussian param
%
% A=max(data);
% x0=mean(xdata);
% s=0.5*(max(xdata)-min(xdata));
% f=A*exp(-((x0)/s).^2);
% plot(xdata,f)
%
%
%
%
%
%
%
%
%
%
% err=inf;
%
% while err>tol
%
%     f=A*exp(-((cr-x0)/s).^2);
%
%     d=cr-f;
%
%     %construct matriz Z
%
%     z1=f/A;
%     z2=2*f.*(pha-x0)/s^2;
%     z3=2*f.*(pha-x0).^2/s^3;
%     z=[z1 z2 z3];
%
%     %update coeficient
%     da=(z.'*z)\(z.'*d);
%     A=A+da(1);
%     x0=x0+da(2);
%     s=s+da(3);
%
%     %calculate error
%
%     err=sum(abs(da./[A;x0;s]));
% end
%
% figure()
% x=linspace(0,4,1000);
% f=A*exp(-((x-x0)/s).^2);
% plot(x,f)
% hold
% plot(pha,cr,'xr')
%
%
%
% %
% %
% % x = cr;
% % y = normpdf(x,0,1);
% %
% % plot(x,y)
% %
% %
% %
% % clearvars -except cr
% % pd_kernel = fitdist(cr,'Kernel')
% % pd_normal = fitdist(cr,'Normal')
% % % x=1:12;
% %
% %  x = min(cr):0.1:max(cr);
% % pdf_kernel = pdf(pd_kernel,x);
% % pdf_normal = pdf(pd_normal,x);
% %
% % plot(x,pdf_kernel,'Color','b','LineWidth',2);
% % hold on;
% % plot(x,pdf_normal,'Color','r','LineStyle',':','LineWidth',2);
% % legend('Kernel Distribution','Normal Distribution','Location','SouthEast');
