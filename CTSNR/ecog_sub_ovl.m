% function [fig]=ecog_sub_ovl(coh_filts,name)

b_in=cell(size(coh_filts,1),1);
b_dur=cell(size(coh_filts,1),1);
for pr=1:size(coh_filts,1)
     pre_b_in=zeros(size(coh_filts{pr,1},1),size(coh_filts{pr,1},2));
    for ct=1:size(coh_filts{pr,1},1)
        sig=coh_filts{pr,1}(ct,:);
        env=abs(hilbert(sig));

%         [onset,offset]=bursts_aligned(env,sig);
%         onset1=onset; clear onset
%         offset1=offset; clear offset
%         
                 [onset1,offset1]=bursts(env);
        onset1=horzcat(onset1{:});
        offset1=horzcat(offset1{:});
            for jj=1:length(onset1)
                pre_b_in(ct,onset1(jj):offset1(jj))=1;
                pre_b_dur(ct,jj)=offset1(jj)-onset1(jj);
            end      
        clear sig env onset1 offset1 
    end
    b_in{pr,1}=pre_b_in;
    b_dur{pr,1}=sum(pre_b_dur,2);
    
    clear pre_b_in pre_b_dur
end

r=0;
perct=[];
prct_all=[];
for iii=1:size(b_in,1)
        probe=squeeze(b_in{iii,1});
  for ct=2:size(probe,1)
      r=r+1;
      ref_b=find(probe(1,:)==1);
      dum=[probe(1,:);probe(ct,:)];
      combi=1;
       [number_ovl,dur_ovl]=overlap(dum,combi); 
       [dur_chance_ovl]=chance_ovl(dum,ref_b,combi);
      
      
      ovlap(r,:)=[(length(ref_b)-dur_ovl) dur_ovl]./length(ref_b).*100;
      ch_ovl(r,:)=dur_chance_ovl;
      
      
      dur_ovlap(r,:)=dur_ovl;
      dur_s(r,:)=[sum(probe(1,:)) sum(probe(ct,:))];
      n_ovlap(1,r)=number_ovl;
      dumi=sum(dum);
      join(r,:)=dumi;
      prct_all(r,:)=(hist(dumi,0:2))./1000;
      clear  dum ref_b
  end
    clear test pre probe_prc probe dur_ovl number_ovl
end


load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','squash','blood','sky','aegean','grey');
if name=='bz'
    color_b={grey,blood};
else
    color_b={grey,aegean};
end

Z=[mean(dur_s) mean(dur_ovlap)]./1000;

fig=figure;
venn(Z,'FaceColor',color_b,'FaceAlpha',{0.6,0.6},'EdgeColor','none')
if name=='bz'
legend({'Ecog'; 'BZ'});
else
legend({'Ecog'; 'SNr'});
end
legend('boxoff')

r=round([(mean(prct_all(:,1))) ((mean(dur_s)-mean(dur_ovlap))./1000) (mean(dur_ovlap)./1000)],1);

str = [num2str(r(1)),'%'];
t = text(-3,2,str);
s = t.FontSize;
t.FontSize = 14;

str = [num2str(r(2)),'%'];
t = text(-0.5,0,str);
s = t.FontSize;
t.FontSize = 14;

str = [num2str(r(3)),'%'];
t = text(3.5,0,str);
s = t.FontSize;
t.FontSize = 14;

str = [num2str(r(4)),'%'];
t = text(1.5,0,str);
s = t.FontSize;
t.FontSize = 14;

set(gca,'visible','off')
fig.Units = 'centimeters';
fig.OuterPosition= [10, 10, 16, 12];



% bar(perct','EdgeColor','none')
% xticklabels({'NO B','NO OVL','OVL'})
% ylabel('accumulated time(s)')
% box('off')

% fig=figure
% bar(mean(prct_all),'EdgeColor','none','FaceColor',color_b{1,2},'FaceAlpha',0.5)
% hold on
% plot([1:3],prct_all,'.')
% % plot(prct_all','.-','LineWidth',2','MarkerSize',15)
% xticklabels({'NO B','NO OVL','OVL'})
% ylabel('Total time (%)')
% box('off')
% set(gca,'FontSize',12)
% fig.Units='centimeters';
% fig.OuterPosition= [10, 10, 12, 10];
% set(fig,'color','w');

fig=figure
bar(mean(ovlap),'EdgeColor','none','FaceColor',color_b{1,2},'FaceAlpha',0.5)
hold on
plot([1:2],ovlap,'.','Color',[0.5 0.5 0.5])
plot([1:2],mean(ch_ovl),'kd','MarkerSize',5,'LineWidth',1)
if (ttest(ovlap(:,1),ch_ovl(:,1)))==1
    plot(1,95,'k*','LineWidth',1)
end
if (ttest(ovlap(:,2),ch_ovl(:,2)))==1
    plot(2,95,'k*','LineWidth',1)
end


% plot((prct_ecog.*100)','.-','LineWidth',2','MarkerSize',15)
xticklabels({'NO OVL','OVL'})
ylabel('Total time EcogB (%)')
box('off')
set(gca,'FontSize',12)
fig.Units='centimeters';
fig.OuterPosition= [10, 10, 12, 10];
set(fig,'color','w');

[mean(ch_ovl) std(ch_ovl)]
[mean(ovlap) std(ovlap)]
(ttest(ovlap(:,1),ch_ovl(:,1)))
(ttest(ovlap(:,2),ch_ovl(:,2)))

% fig=figure
% boxplot(n_ovlap)
% ylim([0 50])
% ylabel('number of overlapping instances')
% if name=='bz'
%     xticklabels({'BZ'});
% else
%     xticklabels({'SNr'});
% end
% box('off')
% set(gca,'FontSize',12)
% fig.Units='centimeters';
% fig.OuterPosition= [10, 10, 5, 10];
% set(fig,'color','w');
% 
%  
% % end