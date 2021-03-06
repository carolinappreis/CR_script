
clear all
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/cleaned_rc12.mat')
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/newnonstim2.mat')

load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\cleaned_rc12_noaddon.mat')
load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\newnonstim10.mat')



main=[1 1 3 1 3 3 3 3 1 1];
%  main=[1 1 1 1 1 1 1 1 1 1];
for pp=1:10
a.ns(pp,:)=nostimout(pp,main(pp),:); a.s(pp,:)=nanmedian(tt1{pp,1}); 
end

% % load ('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\freq_FRC.mat');load ('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\freq_NS.mat')
% load ('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/freq_FRC.mat'); load ('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/freq_NS.mat');
% 
% f.ns=NS(off,:); f.s=ttall(off,:); clearvars -except a f

ref=a.s; %%% max amplitude change vs. max frequecy change
iii=1; %%%%% amp (=0) vs. supressive effect

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
    sub=ref(idma,:);  % iii=0 amplifying effect;
    sub1=a.ns(idma,:);  % iii=0 amplifying effect;
else
    sub=ref(idmi,:);  % iii~=0 supressive effect;
    sub1=a.ns(idmi,:);  % iii=0 amplifying effect;
end

for i=1:size(sub,1);
    
    if iii==0
        phase_peak(1,i)=find(sub(i,:)==max(sub(i,:)));
        p_p_s(1,i)=find(sub1(i,:)==max(sub1(i,:)));
    else
        phase_peak(1,i)=find(sub(i,:)==min(sub(i,:)));
        p_p_s(1,i)=find(sub1(i,:)==min(sub1(i,:)));
    end
    
    if phase_peak(i)==1;
        a_s_al(i,:)=a.s(i,:);
        a_ns_al(i,:)=a.ns(i,:);
        
        %f_s_al(i,:)=f.s(i,:);
        %f_ns_al(i,:)=f.ns(i,:);
        
    else
        a_s_al(i,:)=[a.s(i,phase_peak(i):end) a.s(i,1:phase_peak(i)-1)];
        a_ns_al(i,:)=[a.ns(i,p_p_s(i):end) a.ns(i,1:p_p_s(i)-1)];
        
        %f_s_al(i,:)=[f.s(i,phase_peak(i):end) f.s(i,1:phase_peak(i)-1)];
        %f_ns_al(i,:)=[f.ns(i,phase_peak(i):end) f.ns(i,1:phase_peak(i)-1)];
        
        
        check=[phase_peak(i):size(a_s_al,2) 1:phase_peak(i)-1];
    end
    
    if (phase_peak(i)+5)<= size(a_s_al,2)
        
        a_s2(i,:)= [a.s(i,phase_peak(i)+5:end) a.s(i,1:phase_peak(i)+5-1)];
        a_ns2(i,:)= [a.ns(i,phase_peak(i)+5:end) a.ns(i,1:phase_peak(i)+5-1)];
        %f_s2(i,:)= [f.s(i,phase_peak(i)+5:end) f.s(i,1:phase_peak(i)+5-1)];
        %f_ns2(i,:)= [f.ns(i,phase_peak(i)+5:end) f.ns(i,1:phase_peak(i)+5-1)];
        
    else
        
        dum=(phase_peak(i)+5)-size(a_s_al,2);
        
        a_s2(i,:)= [a.s(i,dum+5:end) a.s(i,1:dum+5-1)];
        a_ns2(i,:)= [a.ns(i,dum+5:end) a.ns(i,1:dum+5-1)];
        %f_s2(i,:)= [f.s(i,dum+5:end) f.s(i,1:dum+5-1)];
        %f_ns2(i,:)= [f.ns(i,dum+5:end) f.ns(i,1:dum+5-1)];
        clear dum
        
    end
    
end
% 
% 
% st=NaN(1,12);
% clear A; A=a_s_al; %b1{f,1};
% clear B; B=a_ns_al; %s1{f,1}(1:size(A,1),:);
% hayriye_c; st(1,:)=stats.prob; st2(1,:)=stats.posclusterslabelmat;
% beg=find(st(1,:)<0.05 & st2(1,:)~=0);
% if ~isempty(beg)
%     sig_rise_all=[beg(1) beg(end)];
%     
% end
% 



%%%%%%%%%%%%%%%%%%%check alignments
% for i=1:10
%     figure(1)
% subplot(1,10,i)
% bar(a.s(i,:))
% 
% 
% % figure(2)
% % subplot(1,10,i)
% % bar(a.ns(i,:))
% 
% figure(3)
% subplot(1,10,i)
% bar(a_s_al(i,:))
% 
% 
% % figure(4)
% % subplot(1,10,i)
% % bar(a_ns_al(i,:))
% end




if kstest(a_s_al(:,1)-a_ns_al(:,1))==1
    
    for i=1:12
        [p1(1,i),h1(1,i)]=signrank(a_s_al(:,i),a_ns_al(:,i));
        a.s_ns=[h1(1) p1(1)];
        
        [p2(1,i),h2(1,i)]=signrank(a_s_al(:,i),a_s2(:,i));
        a.s_180=[h2(1) p2(1)];
        
        [p3(1,i),h3(1,i)]=signrank(a_ns_al(:,i),a_ns2(:,i));
        a.ns_180=[h3(1) p3(1)];
        
        
        
    end
    test_a='wilcoxon';
    
    
else
    
    [p1,h1]=ttest(a_s_al,a_ns_al);
    a.s_ns=[p1(1) h1(1)];
    
    [p2,h2]=ttest(a_s_al,a_s2);
    a.s_180=[p2(1) h2(1)];
    
    [p3,h3]=ttest(a_ns_al,a_ns2);
    a.ns_180=[p3(1) h3(1)];
    test_a='ttest';
end


% clearvars -except a

a.s_ns 

a.s_180

a.ns_180


% for i=1:size(sub,1);
%         new_as_al(i,:)=zscore(a_s_al(i,:));
% end
% new_as_al1=[new_as_al(:,8:12) new_as_al(:,1:7)];
% plot(new_as_al1','-d')
% xlim([1 12])
% xticklabels({'-4','-2','0','2','4',''})
% box('off')
% hold on
% plot(mean(new_as_al1),'k','LineWidth',2)

% clearvars -except test_a test_f f a 