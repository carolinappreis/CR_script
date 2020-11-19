function [fig]=ecog_sub_ovl(coh_filts,name)

b_in=cell(size(coh_filts,1),1);
b_dur=cell(size(coh_filts,1),1);
for pr=1:size(coh_filts,1)
     pre_b_in=zeros(size(coh_filts{pr,1},1),size(coh_filts{pr,1},2));
    for ct=1:size(coh_filts{pr,1},1)
        sig=coh_filts{pr,1}(ct,:);
        env=abs(hilbert(sig));
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
      ref_b=find(probe(ct,:)==1);
      dum=sum([probe(1,:);probe(ct,:)]);
      
      intersect(find(probe(ct,:)==1),find(probe(1,:)==1))
      [number_ovl,dur_ovl]=overlap(dum); 
      ovlap(r,:)=[(length(ref_b)-dur_ovl) dur_ovl]./length(ref_b).*100;
      n_ovlap(1,r)=number_ovl;
      join(r,:)=dum;
      prct_all(r,:)=(hist(dum,0:2))./1000;
      clear  dum ref_b
  end
    clear test pre probe_prc probe dur_ovl number_ovl
end

load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','squash','blood','sky','aegean');
if name=='bz'
    color_b={squash blood};
else
    color_b={sky aegean};
end

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
plot([1:2],ovlap,'.')
% plot((prct_ecog.*100)','.-','LineWidth',2','MarkerSize',15)
xticklabels({'NO OVL','OVL'})
ylabel('Total time EcogB (%)')
box('off')
set(gca,'FontSize',12)
fig.Units='centimeters';
fig.OuterPosition= [10, 10, 12, 10];
set(fig,'color','w');

fig=figure
boxplot(n_ovlap)
ylim([0 50])
ylabel('number of overlapping instances')
if name=='bz'
    xticklabels=({'BZ'});
else
    xticklabels=({'SNr'});
end
box('off')
set(gca,'FontSize',12)
fig.Units='centimeters';
fig.OuterPosition= [10, 10, 5, 10];
set(fig,'color','w');

 
end