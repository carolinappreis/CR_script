%change load file to the region of intrest to check if there is any
%significant different between the change in beta evoked by bursts vs.
%change in beta evoked by surrogates.
% results: only SNr, both short and long bursts show significant windows
clear all
cd ('\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A1_Thal\mat')
load ('SNr_nc_fig','all_contacts','surr2','surr_s')
out3_real= all_contacts;
clear all_contacts
out3_sur= surr2;
clear surr2
color_BZ=[0.5 0 0.5];
color_SNR=[0 0 0.5];
color_surr=[0.3 0.3 0.3];
color_region=color_BZ;

cd ('C:\Users\creis\Documents\GitHub\CRcode\codes_thal\A1_Thal\code')

st=NaN(2,551);
for ind=1:2; clear A; A(1:size(out3_real,2),1:551)=out3_real(ind,:,300:850);
    clear B; B(1:size(out3_sur,1),1:551)=out3_sur(:,300:850);
    hayriye_c; st(ind,:)=stats.prob; st2(ind,:)=stats.posclusterslabelmat;
%     figure(ind)
%     plot(median(B))
%     hold on
%     plot(median(A))
%     plot(st(ind,:),'r.')
%     plot(st2(ind,:),'b.')
    
end

beg=find(st(1,:)<0.05 & st2(1,:)~=0);
if ~isempty(beg)
    beg(1)
    beg(find(diff(beg)>1))
    beg(find(diff(beg)>1)+1)
    beg(end)
    figure(3)
    ax1=subplot(2,2,1)
    patch([beg(1) beg(end) beg(end) beg(1)]-200,[0.2 0.2 0.225 0.225],color_region,'EdgeColor','none')
end
hold on
clear beg
beg=find(st(2,:)<0.05 & st2(2,:)~=0);
if ~isempty(beg)
    beg(1)
    beg(find(diff(beg)>1))
    beg(find(diff(beg)>1)+1)
    beg(end)
    figure(3)
    ax2=subplot(2,2,3)
    patch([beg(1) beg(end) beg(end) beg(1)]-200,[0.2 0.2 0.225 0.225],color_region,'EdgeColor','none')
end
hold on
load ('SNr_nc_fig','med_thal_s')
figure(3)
ax1=subplot(2,2,1)
x=[-500:500];
y1=med_thal_s{1}(1,:);
y2=med_thal_s{1}(2,:);
y3=med_thal_s{1}(3,:);
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_region);
patch([x fliplr(x)], [y1 fliplr(y2)],color_region,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_region,'FaceAlpha',0.2,'EdgeColor','none')
hold on
y1=surr_s(1,:);
y2=surr_s(2,:);
y3=surr_s(3,:);
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_surr);
patch([x fliplr(x)], [y1 fliplr(y2)],color_surr,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_surr,'FaceAlpha',0.2,'EdgeColor','none')
ax2=subplot(2,2,3)
y1=med_thal_s{2}(1,:);
y2=med_thal_s{2}(2,:);
y3=med_thal_s{2}(3,:);
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_region);
patch([x fliplr(x)], [y1 fliplr(y2)],color_region,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_region,'FaceAlpha',0.2,'EdgeColor','none')
hold on
y1=surr_s(1,:);
y2=surr_s(2,:);
y3=surr_s(3,:);
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_surr);
patch([x fliplr(x)], [y1 fliplr(y2)],color_surr,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_surr,'FaceAlpha',0.2,'EdgeColor','none')


% for i =1:79
% plot(A(i,:))
% end