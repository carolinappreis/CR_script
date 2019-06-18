cd('C:\Users\creis\Documents\GitHub\CR_script\Tremor')

cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stimload stim')
load stim
load nostim
d=nanmedian(stimout);

upperthreshold=prctile(nostimout,99.7917); %bonferroni corrected for 12 comparisons
lowerthreshold=prctile(nostimout,0.2083);

find(d>=upperthreshold | d<=lowerthreshold)
bar(d) 
hold on 
rr(1:12)=upperthreshold;
rr2(1:12)=lowerthreshold;
plot(rr,'r--')
plot(rr2,'r--')
xlim([0.5 12.5])
box('off')
xticklabels({'0' '30' '60' '90' '120' '150' '180' '210' '240' '270' '300' '330'})

cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stimload stim')
load ('pha_suffle')
for s=1:size(tt,2)
    for i =1:100;
        yy1=xx(randperm(size(xx,2)) );
        tt2(1:sum(yy1==s),1)=tremor_or2(find(yy1==s));
        tt3(i,s)=nanmedian(tt2,1);
        clear tt2 
    end
end
lg=0:20;
for i=1:12
    tt(1:sum(xx==i),i)=tremor_or2(find(xx==i));
end

rr(1:size(tt,2))=mean(prctile(tt3,95));
rr1(1:size(tt,2))=mean(prctile(tt3,25));
bar(nanmedian(tt))
hold on
plot(rr,'k--','LineWidth',1.5)
plot(rr1,'k--','LineWidth',1.5)
box('off')
title 'Significant phasic-stimulation effect' 
xticklabels({'0' '30' '60' '90' '120' '150' '180' '210' '240' '270' '300' '330'})




figure()
likhood_amp=sum(tt>prctile(tt3,95)| tt>0)./sum(~isnan(tt));
likhood_sup=sum(tt<prctile(tt3,25)| tt<0)./sum(~isnan(tt));
likhood=[likhood_sup ; likhood_amp];
bar(likhood')
title ('Likelihood of significant amplification/supression effect')
xticklabels({'0' '30' '60' '90' '120' '150' '180' '210' '240' '270' '300' '330'})
legend('supression','amplification')
legend('boxoff')
box('off')
