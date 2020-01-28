clear all


iii=0; %%%%% amp (=0) vs. supressive effect
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
% plot(ns_filt_mainax{i,1})
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
        [Pxx_s,F]=pwelch(sup,samplerate,[],samplerate,samplerate);
    else
        Pxx_s=NaN(1,501);
    end
    
    if ~isempty(amp)
        [Pxx_a,F]=pwelch(amp,samplerate,[],samplerate,samplerate);
    else
        Pxx_a=NaN(1,501);
    end
    
    [Pxx_ns,F]=pwelch(ns_filt_mainax{y,1},samplerate,[],samplerate,samplerate);
    
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
    
    
    clear Pxx_nstim
    clear Pxx_a
    clear Pxx_s
end
samplerate=1000;
idx=find(~isnan(Pxx_s_all(:,1)));
S=max([Pxx_s_all(idx,:)]')
A=max([Pxx_a_all(idx,:)]')
NS=max([Pxx_ns_all(idx,:)]')
all=[S ; A ; NS];

figure;
bar(1:3,[median(S) median(A) median(NS)],'FaceColor',[0 0.5 0.5],'EdgeColor','none')
hold on
plot([S; A; NS],'.','MarkerSize',5,'Color','k')
xlim([0 4])
xticklabels({'AMP','SUP','NS'})
ylim([0 31])
yticks([1:2:30])
breakyaxis([8 28]);
box('off')
ylabel('Peak pwelch')
%  BreakPlot([median(S) median(A) median(NS)],[1:13,27:30],13,27,'Patch');

% figure
% plot(all,'k')
% hold on
% plot([S; A; NS],'.','MarkerSize',5,'Color','k')
% xlim([0 4])
% xticks([13])
% xticklabels({'AMP','SUP','NS'})
% ylim([0 31])
% yticks([1:2:30])
% box('off')



[p,h]=ttest(A,S)
[p,h]=ttest(A,NS)
[p,h]=ttest(S,NS)


