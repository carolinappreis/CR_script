clear all


iii=0; %%%%% amp (=0) vs. supressive effect
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/seg_env_perphase.mat')
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


sup=[];
amp=[];

for y=1:size(Sa,1)
    for yy=1:size(Sa,2)
        dum=seg_filt{y,yy};
        dum2=reshape(dum',1,(size(dum,1)*size(dum,2)));
        if Sa(y,yy)>0
            amp=[amp dum2];
        else
            sup=[sup dum2];
        end
        clear dum dum2
    end
end
samplerate=1000;
[Pxx_s,F]=pwelch(sup,samplerate,[],samplerate,samplerate);
[Pxx_a,F]=pwelch(amp,samplerate,[],samplerate,samplerate);

plot(F,Pxx_s,'Color','b','LineWidth',2)
hold on
plot(F,Pxx_a,'Color','r','LineWidth',2)
xlim([0 15])
legend ('sup','amp')
legend('boxoff')
box('off')

[p,h]=ttest2(amp,sup)
