function [fig]=subctx_ovl(coh_filts,name)

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
for iii=1:size(b_in,1)
    if size(b_in{iii,1},1)>=3
        r=r+1;
        probe=squeeze(b_in{iii,1}(2:end,:));
        probe1=sum(probe);
        ovl=hist(probe1,0:max(probe1));
        perct(r,:)=[ovl(1)./1000 ovl(2)./1000 sum(ovl(3:end))./1000];
        clear ovl probe bins s_test
    end
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

fig=figure;
plot(1:3,perct','.-','LineWidth',2','MarkerSize',15)
hold on
bar(1:3,median(perct),'EdgeColor','none','FaceColor',color_b{1,2},'FaceAlpha',0.5);
xlim([0 4])
xticks([1:3])
xticklabels({'NO B','NO OVL','OVL'})
ylabel('Total time (%)')
box('off')
set(gca,'FontSize',12)
fig.Units='centimeters';
fig.OuterPosition= [10, 10, 12, 10];
set(fig,'color','w');

end