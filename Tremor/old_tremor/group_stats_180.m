clear all
load ('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\F_group.mat')
% load ('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/F_group.mat')

LS1=squeeze(mean(LS,2));

for i=1:size(S,1);
    phase_peak(1,i)=find(S(i,:)==max(S(i,:)));
    
    if phase_peak(i)==1;
        s1(i,:)=S(i,:);
        ns1(i,:)=NS(i,:);
       
    else
        s1(i,:)=[S(i,phase_peak(i):end) S(i,1:phase_peak(i)-1)];
        ns1(i,:)=[NS(i,phase_peak(i):end) NS(i,1:phase_peak(i)-1)];
    end
    
    
    if (phase_peak(i)+5)<= size(S,2)
        
    s2(i,:)= [S(i,phase_peak(i)+5:end) S(i,1:phase_peak(i)+5-1)];
    ns2(i,:)= [NS(i,phase_peak(i)+5:end) NS(i,1:phase_peak(i)+5-1)];
    
    else
        
    dum=(phase_peak(i)+5)-size(S,2);
        
    s2(i,:)= [S(i,dum+5:end) S(i,1:dum+5-1)];
    ns2(i,:)= [NS(i,dum+5:end) NS(i,1:dum+5-1)];
    clear dum
            
    end

end
[p,h1]=ttest(s1,s2);
s_op_phase=([find(h1<(0.05/size(s1,1))) h1(find(h1<(0.05/size(s1,1)))) h1(1)])
[p,h2]=ttest(ns1,ns2);
ns_op_phase=([find(h2<(0.05/size(s1,1))) h2(find(h2<(0.05/size(s1,1)))) h2(1)])

