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

figure()
for i=1:10
subplot(5,2,i)
plot(seg_filt{i,1})
% plot(ns_filt_mainax{i,1})
end

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

f1=figure(1);
subplot(1,2,1)
bar(1:3,[median(S) median(A) median(NS)],'FaceColor',[0.13,0.22,0.27],'FaceAlpha',0.4,'EdgeColor','none')
hold on
plot([S; A; NS],'.','MarkerSize',10,'Color','k')
yticks([0:0.4:2])
xlim([0 4])
xticklabels({'AMP','SUP','NS'})
box('off')
ylabel('Theta power peak')
% f1.Units = 'centimeters';
% f1.OuterPosition= [10, 10, 8, 8];
set(gca,'FontSize',12)
set(f1,'color','w');
f1=figure(1);
subplot(1,2,2)
for i=1:8
cl= rand(1,3);
plot(all(:,i),'Color',cl,'LineWidth',1)
hold on
plot([S(i); A(i); NS(i)],'.','MarkerSize',15,'Color',cl)
xlim([0 4])
xticks([1:3])
xticklabels({'AMP','SUP','NS'})
box('off')
ylabel('Theta power peak')
yticks([0:0.2:1.6])
end
set(gca,'FontSize',12)
set(f1,'color','w');

[p,h]=ttest(A,S)
[p,h]=ttest(A,NS)
[p,h]=ttest(S,NS)


