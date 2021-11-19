 function [f1,y]=trem_prop(out,match_ax,spiral)  %%% gaussian fit to tremor properties & correlation with arc features


load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','stone');cl=stone;

for iii=1:size(out.start_c,1)
    

    
    y=out.amp_n_bins(iii,:);
%     width(iii,:)=sum(~isnan(y));
    y(1,find(isnan(y)))=0;
   
    x=out.bins(1:end-2);
    pre_fwhm=x(find(y>0.5));
    fwhm(iii,:)=pre_fwhm(end)-pre_fwhm(1);
    
    f1=figure(1)
    subplot(1,4,iii)
    bar(x,y,'FaceColor',cl,'EdgeColor',cl)
    hold on
    [rsg,rsg_g,rsg_o]=gauss_fit2(x',y');
     h2=plot(rsg);
%         pd = makedist('normal','mu',rsg.b,'sigma',rsg.c);

    
    
    [fitobj,goodness,output] = WBFit(x', y');
    h2=plot(fitobj);
    pd = makedist('Weibull','a',fitobj.a,'b',fitobj.b);
    var1(iii,:)=iqr(pd);
    var2(iii,:)=std(pd)./mean(pd);

%   [fitobj,goodness,output]=fit(x',y','weibull')
%   h2=plot(fitobj);
    set(h2,'LineWidth',1.5,'Color','r')
    %     [rsg,rsg_g,rsg_o] = raylFit(x, y)
    %    [rsg,rsg_g,rsg_o] = WBFit(x,y);
    %     ylim([0 1.5])
    xticks([1:2:14])
    xticklabels({'2','3','4','5','6','7','8'});
    ylabel({'Absolute amplitude (\muV^2)'})
    xlabel('Frequency (Hz)')
    set(gca,'FontSize',12)
    box('off')
    legend('off')
    title(sprintf('patient %d',(iii)))
%      cv(iii,:)=rsg.c./rsg.b;
%     clear y rsg rsg_o rsg_g
    set(gca,'FontSize',12)
    f1.Units = 'centimeters';
    f1.OuterPosition= [10, 10, 25, 10];
    set(f1,'color','w');
    
    dns=prctile((squeeze(out.mod_amp{iii,1}(match_ax(1,iii,1),:))),0.2083);
    d1=nanmedian(squeeze(out.mod_amp{iii,2}{match_ax(2,iii,1),1}));
    arc.value(iii,1,1)=nanmean(abs(d1))*100;
    arc.value(iii,2,1)=(max(d1)-min(d1));
    if min(d1)<0
        arc.value(iii,3,1)=-(min(d1)*100);
    else
        arc.value(iii,3,1)=NaN;
    end
    arc.value(iii,4,1)=-(min(d1)-dns)
    arc.label={'mean effect';'effect range';'supressive effect'; 'deviation from natural change (%)' };
    clear d1 dns
end


f1=figure(2)
%subplot(1,size(arc.value,2),ii)
ii=4;
dum=find(~isnan((arc.value(:,ii))));
%  x=cv(dum);
% x=fwhm(dum);
%         x=width(dum);
          x=var2(dum);

if spiral==0
    y=arc.value(dum,ii);
else
    load('/Users/Carolina/Documents/GitHub/CR_script/DBS_tremor/DBS_paper/all_cond_spiral.mat','Ppeak')
    y=abs((Ppeak(:,3)-Ppeak(:,1))./-Ppeak(:,1));
    x=x(2:end);
end
[fitresult] = fit( x, y, 'poly1' );
h=plot(fitresult,x,y);
legend('off')
hold on
plot(x,y,'k.','MarkerSize',10,'HandleVisibility','off');
[r,p]=corrcoef(x,y); r=vpa(round(r,2));p=vpa(round(p,2));
title(sprintf('r = %s',r),'FontSize',12)
% 
% text(0.16,0.25,sprintf('r = %s',r),'FontSize',12)
% text(0.16,0.20,sprintf('p = %s',p),'FontSize',12)

box('off')

 xlabel('IQ')
% xlabel('fwhm')
%  xlabel('coefficient of variation')
% ylabel(sprintf( arc.label{ii,1}))
ylabel('tremor suppresion')

set(gca,'FontSize',12)
f1.Units ='centimeters';
f1.OuterPosition= [10, 10, 8, 10];
set(f1,'color','w');


 end