clear all
load ('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\A_group.mat')
% load ('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/A_group.mat')
a.ls=LS; a.ns=NS; a.s=S; clearvars -except a
load ('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\F_group.mat')
% load ('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/F_group.mat')
f.ls=LS; f.ns=NS; f.s=S; clearvars -except a f
a.ls1=squeeze(mean(a.ls,2));
f.ls1=squeeze(mean(f.ls,2));

ref=a.s; %%% max amplitude change vs. max frequecy change
iii=0; %%%%% amp (=0) vs. supressive effect

idmi=[];
idma=[];
for i =1:size(ref,1)
    if (numel(find(ref(i,:)<0)))~=0 % all -except all amplifying subjects
        idmi=[idmi i];
    end
    if (numel(find(ref(i,:)>0)))~=0 % all -except all supressive subjects
        idma=[idma i];
    end
end


if iii==0;
    sub=ref(idma,:)  % iii=0 amplifying effect;
else
    sub=ref(idmi,:)  % iii~=0 supressive effect;
end


for i=1:size(sub,1);
    
    if iii==0
        phase_peak(1,i)=find(sub(i,:)==max(sub(i,:)));
    else
        phase_peak(1,i)=find(sub(i,:)==min(sub(i,:)));
    end
    
    if phase_peak(i)==1;
        a_s_al(i,:)=a.s(i,:);
        a_ns_al(i,:)=a.ns(i,:);
        a_ls_al(i,:)=a.ls1(i,:);
        
        f_s_al(i,:)=f.s(i,:);
        f_ns_al(i,:)=f.ns(i,:);
        f_ls_al(i,:)=f.ls1(i,:);
    else
        a_s_al(i,:)=[a.s(i,phase_peak(i):end) a.s(i,1:phase_peak(i)-1)];
        a_ns_al(i,:)=[a.ns(i,phase_peak(i):end) a.ns(i,1:phase_peak(i)-1)];
        a_ls_al(i,:)= [a.ls1(i,phase_peak(i):end) a.ls1(i,1:phase_peak(i)-1)];
        
        
        f_s_al(i,:)=[f.s(i,phase_peak(i):end) f.s(i,1:phase_peak(i)-1)];
        f_ns_al(i,:)=[f.ns(i,phase_peak(i):end) f.ns(i,1:phase_peak(i)-1)];
        f_ls_al(i,:)= [f.ls1(i,phase_peak(i):end) f.ls1(i,1:phase_peak(i)-1)];
        
        check=[phase_peak(i):size(a_s_al,2) 1:phase_peak(i)-1];
    end
    
        if (phase_peak(i)+5)<= size(a_s_al,2)
        
    a_s2(i,:)= [a.s(i,phase_peak(i)+5:end) a.s(i,1:phase_peak(i)+5-1)];
    a_ns2(i,:)= [a.ns(i,phase_peak(i)+5:end) a.ns(i,1:phase_peak(i)+5-1)];
    f_s2(i,:)= [f.s(i,phase_peak(i)+5:end) f.s(i,1:phase_peak(i)+5-1)];
    f_ns2(i,:)= [f.ns(i,phase_peak(i)+5:end) f.ns(i,1:phase_peak(i)+5-1)];
    
    else
        
    dum=(phase_peak(i)+5)-size(a_s_al,2);
        
    a_s2(i,:)= [a.s(i,dum+5:end) a.s(i,1:dum+5-1)];
    a_ns2(i,:)= [a.ns(i,dum+5:end) a.ns(i,1:dum+5-1)];
    f_s2(i,:)= [f.s(i,dum+5:end) f.s(i,1:dum+5-1)];
    f_ns2(i,:)= [f.ns(i,dum+5:end) f.ns(i,1:dum+5-1)];
    clear dum
            
    end
    
end

tests={'signrank','ttest'};
if kstest(a_s_al(:,1)-a_ns_al(:,1))==1
    stat=tests{1};
else
    stat=tests{2};
end
%%%% changeeeeeeeeeeeeeeeee
disp(eval(stat(a_s_al(:,1),a_ns_al(:,1))))

[p,h1]=stat(a_s_al,a_ns_al);
a.s_ns=([find(h1<(0.05/size(a_s_al,2)))  h1(find(h1<(0.05/size(a_s_al,2)))) h1(1)])

[p,h2]=ttest(a_s_al,a_ls_al);
a.s_rs=([find(h2<(0.05/size(a_s_al,2))) h2(find(h2<(0.05/size(a_s_al,2)))) h2(1)])

[p,h3]=ttest(a_s_al,a_s2);
a.s_180=([find(h3<(0.05/size(a_s_al,2))) h3(find(h3<(0.05/size(a_s_al,2)))) h3(1)])

[p,h4]=ttest(a_s_al,a_s2);
a.ns_180=([find(h4<(0.05/size(a_s_al,2))) h4(find(h4<(0.05/size(a_s_al,2)))) h4(1)])

clear p h1 h2 h3 h4


[p,h1]=ttest(f_s_al,f_ns_al);
f.s_ns=([find(h1<(0.05/size(f_s_al,2))) h1(find(h1<(0.05/size(a_s_al,2)))) h1(1)])

[p,h2]=ttest(f_s_al,f_ls_al);
f.s_rs=([find(h2<(0.05/size(f_s_al,2)))  h2(find(h2<(0.05/size(a_s_al,2)))) h2(1)])

[p,h3]=ttest(f_s_al,f_s2);
f.s_180=([find(h3<(0.05/size(f_s_al,2))) h3(find(h3<(0.05/size(f_s_al,2)))) h3(1)])

[p,h4]=ttest(f_s_al,f_s2);
f.ns_180=([find(h4<(0.05/size(f_s_al,2))) h4(find(h4<(0.05/size(f_s_al,2)))) h4(1)])

clearvars -except a f
