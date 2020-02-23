clear all
close all
%
% load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\cleaned_rc12_noaddon.mat')
% load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\newnonstim10.mat')
% load('C:\Users\creis\Documents\GitHub\CR_script\colour_pal.mat','blushred','squash');
% cl=blushred;
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/cleaned_rc12_noaddon.mat')
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/newnonstim10.mat')
cl=[0.5 0.5 0.5];

main=[1 1 3 1 3 3 3 3 1 1];
for pp=1:size(tt1,1)
    for kk=1:size(tt1,2)
        
        curve=nanmedian(tt1{pp,kk},1);
        SO=repmat(curve,1,3);
        %         SOUT=repmat((squeeze(nostimout(pp,kk,:)))',1,3);
        for i=size(tt1{pp,kk},2)+1:size(tt1{pp,kk},2)*2
            smooth_c(1,i-12)=sum(SO(1,(i-1:i+1)))./length(SO(1,(i-1:i+1)));
            %             smooth_ns(1,i-12)=sum(SOUT(1,(i-1:i+1)))./length(SOUT(1,(i-1:i+1)));
        end
        
        
        smoo_all(pp,kk,:)=smooth_c;  clear smooth_c
        %         nsmoo_all(pp,kk,:)=smooth_ns; clear smooth_ns
        
    end
    
    smoo_main(pp,:)=squeeze(smoo_all(pp,1,:));
    %     nsmoo_main(pp,:)=squeeze(nsmoo_all(pp,1,:));
    s_main(pp,:)=nanmedian(tt1{pp,1},1);
    %     ns_main(pp,:)=squeeze(nostimout(pp,1,:));
end



aic=[];
for i=1:size(s_main,1);
    y= s_main(i,:);  
    f1=figure(1)
    subplot(2,5,i)
    bar(y,'FaceColor',cl,'EdgeColor',cl,'HandleVisibility','off')
    hold on
    [rsg,rsg_g,rsg_o]=gauss_fit(y);
    N=rsg_o.numobs;
    SS=(rsg_g.rmse)^2  ;
    K=rsg_o.numparam;
    DFE=rsg_g.dfe;
    Var = DFE/N*SS;
%     Ll = -( sum((y'-rsg(1:length(y))).^2)/Var + N*log(2*pi) + N*log(Var))/2;
     Ll = sum((-((y'-rsg(1:length(y))).^2))./(2*(SS^2)));
    aic(3) = real(aicbic((Ll),K));
    clear N SS K Ll 
    
    [rss,rss_g,rss_o]=sin_fit(y);
    N=rss_o.numobs;
    SS=(rss_g.rmse)^2 ;
    K=rss_o.numparam;
    DFE=rss_g.dfe;
    Var = DFE/N*SS;
%     Ll = -( sum((y'-rss(1:length(y))).^2)/Var + N*log(2*pi) + N*log(Var))/2;
    Ll = sum((-((y'-rss(1:length(y))).^2))./(2*(SS^2)));
    aic(2) = real(aicbic((Ll),K));
    clear N SS K Ll 

    
    x=0:length(y)-1;
    mdk= fitlm(x,y,'constant');
    if sum(diff(mdk.Fitted))==0
    yline(mdk.Fitted(1),'k-.','LineWidth',1.5)
    else
    plot(mdk.Fitted,'k-.','LineWidth',1.5)
    end
    
%     aicbic(mdk.LogLikelihood,mdk.NumEstimatedCoefficients)
%     mdk.ModelCriterion.AIC
%     mdk.ModelCriterion.AICc
    N=mdk.NumObservations;
    SS=mdk.MSE;
    K=mdk.NumEstimatedCoefficients;
    DFE=mdk.DFE;
 
    mdk= fitlm(x,y,'constant');
    N=mdk.NumObservations;
    SS=mdk.MSE;
    K=mdk.NumEstimatedCoefficients;
    DFE=mdk.DFE;
    Var = DFE/N*SS;
%     Ll = -( sum((y'-mdk.Fitted).^2)/Var + N*log(2*pi) + N*log(Var))/2;
    Ll = sum((-((y'-mdk.Fitted).^2))./(2*(SS^2)));
    aic(1) = -2*Ll + 2*K;
%     real(aicbic(Ll,K));
    
    win(i,1)=find(aic==(min(aic)));
    clear aic

    xticklabels({'0','30','60','90','120','150','180','210','240','270','300','330'});
    ylabel('Change in tremor severity')
    xlabel('Stimulation phase (degrees)')
    set(gca,'XTickLabelRotation',45)
    set(gca,'FontSize',12)
    box('off')
    legend({'gaussian fit','sinusoidal fit','constant fit'},'Location','northwest')
    legend('boxoff')
    set(legend,'FontSize',9)
legend('off')
    
end

f1.Units = 'centimeters';
f1.OuterPosition= [10, 10, 55, 15];
set(f1,'color','w');


% close all
% for i=1:size(s_main,1)
%
%    y= s_main(i,:);
%    f1=figure(1)
%     subplot(2,5,i)
%     bar(y,'FaceColor',cl,'EdgeColor',cl)
%     hold on
%     poli_fit(y);
%    sin_fit(y);
%    gauss_fit(y);
% legend('off')
%    xlabel('');ylabel('')
%    box('off')
%  title({['p k-sin:' num2str(f_test(1,i))];['p k-gaus:' num2str(f_test(2,i))]})
% end
%

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
