
%trial count
load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_cluster_out.mat');
for iii=1:4
dum=[];dum1=[];dum2=[];
for i=1:12
dum=[ dum sum(~isnan(out.change_c{iii,2}{1,1}(:,i)))];
dum1=[dum1 sum(~isnan(out.arc1{iii,2}{1,1}(:,i)))];
dum2= [dum2 sum(~isnan(out.arc2{iii,2}{1,1}(:,i)))];
end
out.trials(iii,1)= round(mean(dum));
out.trials(iii,2)= round(mean(dum1));
out.trials(iii,3)= round(mean(dum2));
end
