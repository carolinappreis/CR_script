clear all

% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/seg_env_perphase.mat')
load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\NS_filt_mainax.mat')
load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\seg_env_perphase.mat')

cohort=[2 3 4 5 8 10 11 13 16 17];
dum=intersect(iiii,cohort);

pt=[];
for i=1:length(dum)
    pt=[pt find(iiii==dum(i))];
end

for pp=1:size(pt,2)
    Sa(pp,:)=nanmedian(tt1{pt(pp),1})  ;
end

% figure()
% for i=1:10
% subplot(5,2,i)
% plot(seg_filt{i,1})
% % plot(ns_filt_mainax{i,1})
% end

for y=1:size(Sa,1)
    sup=[];
    amp=[];
    for yy=1:size(Sa,2)
        dum=seg_filt{y,yy};
        dum2=reshape(dum',1,(size(dum,1)*size(dum,2)));
        if Sa(y,yy)>0
            amp=[amp dum2];
        else
            sup=[sup dum2];
        end
    end
    
    clear dum dum2
    
    samplerate=1000;
    if ~isempty(sup)
        [Pxx_s,F]=pwelch(sup,samplerate*2,[],[],samplerate);
    else
        Pxx_s=NaN(1,1025);
    end
    
    if ~isempty(amp)
        [Pxx_a,F]=pwelch(amp,samplerate*2,[],[],samplerate);
    else
        Pxx_a=NaN(1,1025);
    end
    
    [Pxx_ns,F]=pwelch(ns_filt_mainax{y,1},samplerate*2,[],[],samplerate);
    
    Pxx_s_all(y,:)=Pxx_s; 
    Pxx_a_all(y,:)=Pxx_a; 
    Pxx_ns_all(y,:)=Pxx_ns; 
    
    
    
%     plot(F,Pxx_s,'Color','b','LineWidth',2)
%     hold on
%     plot(F,Pxx_a,'Color','r','LineWidth',2)
%     plot(F,Pxx_ns,'Color','g','LineWidth',2)
%     xlim([0 15])
%     legend ('sup','amp','NS')
%     legend('boxoff')
%     box('off')
%     
%     
    clear Pxx_nstim
    clear Pxx_a
    clear Pxx_s
end

load('C:\Users\creis\Documents\GitHub\CR_script\colour_pal.mat','blushred','aegean','stone','squash');
samplerate=1000;
idx=find(~isnan(Pxx_s_all(:,1)));
dumm={'Pxx_a_all';'Pxx_s_all';'Pxx_ns_all'};
color_b1={blushred; aegean; stone};

for i=1:size(dumm,1)
    cr=eval(dumm{i,1});
    cr=cr(idx,:);
    ASNS(i,:)=max(cr');
    
    for jo=1:size(cr,1)
        fASNS(i,jo,1)= F(find(cr(jo,:)==(max(cr(jo,:)))));
        seg=[((find(cr(jo,:)==max(cr(jo,:))))-3):((find(cr(jo,:)==max(cr(jo,:))))+3)];
        sASNS(i,jo,:)=cr(jo,seg); clear seg 
    end

    clear cr
end

    [p,h]=ttest(ASNS(1,:),ASNS(2,:))
    [p,h]=ttest(ASNS(1,:),ASNS(3,:))
    [p,h]=ttest(ASNS(2,:),ASNS(3,:))
    [p,h]=ttest(fASNS(1,:),fASNS(2,:))
    [p,h]=ttest(fASNS(1,:),fASNS(3,:))
    [p,h]=ttest(fASNS(2,:),fASNS(3,:))
    
    
close all
f1=figure(1);
color_b1={blushred; aegean; stone};
subplot(1,3,1)
for i=1:3 
    dr=squeeze(sASNS(i,:,:));
y2=median(dr,1);
y1=y2+std(dr,1);
y3=y2-std(dr,1);
time=1:length(y2);
patch([time fliplr(time)], [y1 fliplr(y2)],[color_b1{i,1}],'FaceAlpha',[0.15],'EdgeColor','none','HandleVisibility','off')
patch([time fliplr(time)], [y2 fliplr(y3)],[color_b1{i,1}],'FaceAlpha',[0.15],'EdgeColor','none','HandleVisibility','off')
plot(time,y2,'LineStyle','-', 'LineWidth',1.5,'Color',color_b1{i,1})
set(gca,'FontSize',12)
box('off')
% ylim([0 1.6])
% yticks([0:0.4:2])
xticks([1:1:7])
xticklabels({'-1.5','-1','-0.5','peak','0.5','1','1.5'})
clear y2 
hold on
xlabel('Frequency')
ylabel('Tremor PSD')
hold on
clear y2 y1 y3
end
legend({'AMP','SUP','NS'})
legend('boxoff')
set(gca,'FontSize',12)
legend('boxoff')

subplot(1,3,2)
bar(1:3,[median(fASNS(1,:)) median(fASNS(2,:)) median(fASNS(3,:))],'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.3,'EdgeColor','none')
hold on
subplot(1,3,3)
bar(1:3,[median(ASNS(1,:)) median(ASNS(2,:)) median(ASNS(3,:))],'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.3,'EdgeColor','none')
hold on
for i=1:size(fASNS,2)
cl= rand(1,3);
subplot(1,3,2)
plot(fASNS(:,i),'Color',cl,'LineWidth',1)
hold on
plot([fASNS(1,i); fASNS(2,i); fASNS(3,i)],'.','MarkerSize',15,'Color',cl)
% ylim([0 1.6])
% yticks([0:0.4:2])
xlim([0 4])
xticklabels({'AMP','SUP','NS'})
box('off')
ylabel('Tremor frequency peak')
set(gca,'FontSize',12)

subplot(1,3,3)
plot(ASNS(:,i),'Color',cl,'LineWidth',1)
hold on
plot([ASNS(1,i); ASNS(2,i); ASNS(3,i)],'.','MarkerSize',15,'Color',cl)
% ylim([0 1.6])
% yticks([0:0.4:2])
xlim([0 4])
xticklabels({'AMP','SUP','NS'})
box('off')
ylabel('Tremor power peak')

end
f1.Units ='centimeters';
f1.OuterPosition= [10, 10, 30, 20];
set(f1,'color','w');






