clear all


iii=1; %%%%% amp (=0) vs. supressive effect
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/cleaned_rc12_noaddon.mat')
 load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\cleaned_rc12_noaddon.mat')

cohort=[2 3 4 5 8 10 11 13 16 17];
dum=intersect(iiii,cohort);

pt=[];
for i=1:length(dum)
    pt=[pt find(iiii==dum(i))];
end

for pp=1:size(pt,2)
  Sa(pp,:)=nanmedian(tt1{pt(pp),1})  ;
end


sm=[Sa Sa Sa];
for ii=1:size(sm,1)
    for i=size(Sa,2)+1:size(Sa,2)*2
        a.s(ii,i-12)=sum(sm(ii,(i-1:i+1)))./length(sm(ii,(i-1:i+1)));
    end
end



ref=a.s; %%% max amplitude change vs. max frequecy change

idmi=[];
idma=[];
for i =1:size(ref,1)
    if (numel(find(ref(i,:)<0)))~=0 % all -except all amplifying subjects
        idmi=[idmi i];
    end
    if (numel(find(ref(i,:)>0)))~=0 % all -except all supressive subjects
        idma=[idma i];
    end
end

if iii==0;
    idx=idma;  % iii=0 amplifying effect;
else
    idx=idmi;  % iii~=0 supressive effect;
end

sub=ref(idx,:);

for i=1:length(idx);
    if iii==0
        phase_peak(1,i)=find(sub(i,:)==max(sub(i,:)));
    else
        phase_peak(1,i)=find(sub(i,:)==min(sub(i,:)));
    end
    
    if phase_peak(i)==1;
        a_s_al(i,:)=a.s(idx(i),:);        
        
    else
        a_s_al(i,:)=[a.s(idx(i),phase_peak(i):end) a.s(idx(i),1:phase_peak(i)-1)];        
        check=[phase_peak(i):size(a_s_al,2) 1:phase_peak(i)-1];
    end
    
end


% close all

    metric1=[a_s_al(:,8:12) a_s_al(:,1:7)];
%    load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','squash')
     load('C:\Users\creis\Documents\GitHub\CR_script\colour_pal.mat','blushred','aegean','stone','squash');
    if iii==0
     cl=blushred;
    cl1=squash;
    
    else
       cl=aegean;
    cl1=stone;    
    end
    
% 
% f1=figure(1)
% plot(metric1','Color',cl1)
% ylim([-0.4 0.4])
% xlim([1 12])
% xticks([ 3 6 9 12])
% xticklabels({'-90','0','90','180'})
% box('off')
% hold on
% plot(mean(metric1),'k','LineWidth',4,'Color', cl)
% yline(0,'k', 'LineWidth',1,'LineStyle','--')
% ylabel({'Change in tremor severity'})
% set(gca,'FontSize',14)
% 
% 
% 
% 
% 
% f2=figure(2)
% bar(mean(metric1),'LineWidth',1,'FaceColor',cl,'EdgeColor',cl)
% ylim([-0.4 0.4])
% xlim([0 13])
% xticks([ 3 6 9 12])
% xticklabels({'-90','0','90','180'})
% box('off')
% box('off')
% ylabel({'Change in tremor severity'})
% set(gca,'FontSize',14)
% 
% 
% if iii==0
%  
%     
%     f1=figure(1)
%     xlabel({'Alignment to phase with max';'tremor amplification'})
%     f2=figure(2)
%     xlabel({'Alignment to phase with max';'tremor amplification'})
%     
% else
%     
%     f1=figure(1)
%     xlabel({'Alignment to phase with max';'tremor supression'})
%     f2=figure(2)
%     xlabel({'Alignment to phase with max';'tremor supression'})
%   
% end
% 
% f1.Units = 'centimeters';
% f1.OuterPosition= [10, 10, 12, 12];
% box('off')
% set(f1,'color','w');
% 
% f2.Units = 'centimeters';
% f2.OuterPosition= [10, 10, 12, 12];
% box('off')
% set(f2,'color','w');
% 
figure(10)
subplot(1,2,2)
plot(metric1','Color',cl1)
ylim([-0.4 0.4])
xlim([1 12])
xticks([ 3 6 9 12])
xticklabels({'-90','0','90','180'})
box('off')
hold on
plot(mean(metric1),'k','LineWidth',4,'Color', cl)
yline(0,'k', 'LineWidth',1,'LineStyle','--')
ylabel({'Change in tremor severity'})
xlabel({'Alignment to phase with max';'tremor supression'})
set(gca,'FontSize',14)


% cr=mean(metric1);
% clearvars -except cr
% cr1=mean(metric1);
% close all
% y1=plot(cr,cr1,'k*')
% box('off')
% c2=corrcoef(cr',cr1')
% y3=lsline;
% legend(y3,[num2str(c2)],'box','off')