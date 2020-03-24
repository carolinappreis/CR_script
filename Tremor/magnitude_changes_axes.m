clear all
close all
%  load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\cleaned_rc12_noaddon.mat')
 load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/cleaned_rc12_noaddon.mat')

for i=1:size(tt1,1)
    for ax =1:size(tt1,2)

   abs_change(i,ax,:)=(abs(nanmedian(tt1{i,ax})));
   sum_change(i,ax)=sum(abs(nanmedian(tt1{i,ax})));
    end
end
   


for i=1:size(abs_change)
% [p] = anova1((squeeze(abs_change(i,:,:))'));

    [p,h]=signrank((squeeze(abs_change(i,2,:))'),(squeeze(abs_change(i,1,:))'),0.05./2); 
 

    [p1, h1]=signrank((squeeze(abs_change(i,3,:))'),(squeeze(abs_change(i,1,:))'),0.05./2); 

close all
hr(i,:)=[h h1];; clear n n1
pr(i,:)=[p p1];
end


non_phase_sum(1:10,1:3)=sum(abs_change,3);
bar(non_phase_sum)
box('off')
xlabel('Subjects')
ylabel('Magnitude of modulation') 


for i=1:10
    for aux=2:3
        [h,p]=ttest(squeeze(abs_change(i,1,:)),squeeze(abs_change(i,aux,:)))
       pp(i,aux-1,1)=p;
       answ(i,aux-1,1)=p<0.025;
    end
end