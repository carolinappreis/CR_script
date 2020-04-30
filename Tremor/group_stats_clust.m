

clear all
clear all; close all;
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/data_matchcluster')
load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','squash','aegean','stone');

ref=squeeze((tt1(:,1,:)));
ref_ns=squeeze((nostimout(:,1,:)));

for cr=2
    m=[0 1];
    clearvars -except ref ref_ns cr m f1 blushred squash aegean stone res check
    iii=m(cr); %%%%% amp (=0) vs. supressive effect
    
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
    else
        sub=ref(idmi,:);  % iii~=0 supressive effect;
    end
    
    for i=1:size(sub,1);
        
        if iii==0
            phase_peak(1,i)=find(sub(i,:)==max(sub(i,:)));
        else
            phase_peak(1,i)=find(sub(i,:)==min(sub(i,:)));
        end
        
        if phase_peak(i)==1;
            a_s_al(i,:)=ref(i,:);
            a_ns_al(i,:)=ref_ns(i,:);
            
        else
            a_s_al(i,:)=[ref(i,phase_peak(i):end) ref(i,1:phase_peak(i)-1)];
            a_ns_al(i,:)=[ref_ns(i,phase_peak(i):end) ref_ns(i,1:phase_peak(i)-1)];
            
            
            check(i,:)=[phase_peak(i):size(a_s_al,2) 1:phase_peak(i)-1];
        end
        
            if (phase_peak(i)+5)<= size(a_s_al,2)
        
                a_s2(i,:)= [ref(i,phase_peak(i)+5:end) ref(i,1:phase_peak(i)+5-1)];
                a_ns2(i,:)= [ref_ns(i,phase_peak(i)+5:end) ref_ns(i,1:phase_peak(i)+5-1)];
            else
        
                dum=(phase_peak(i)+5)-size(a_s_al,2);
        
                a_s2(i,:)= [ref(i,dum+5:end) ref(i,1:dum+5-1)];
                a_ns2(i,:)= [ref_ns(i,dum+5:end) ref_ns(i,1:dum+5-1)];
                clear dum
        
            end
        
    end
    
    
    
    if kstest(a_s_al(:,1)-a_ns_al(:,1))==1
        
        for i=1:12
            [p1(1,i),h1(1,i)]=signrank(a_s_al(:,i),a_ns_al(:,i));
            ref_nss=[h1(1) p1(1)];
            
            [p2(1,i),h2(1,i)]=signrank(a_s_al(:,i),a_s2(:,i));
            ref_180=[h2(1) p2(1)];
            
            [p3(1,i),h3(1,i)]=signrank(a_ns_al(:,i),a_ns2(:,i));
            ref_ns_180=[h3(1) p3(1)];
            
            
            
        end
        test_a='wilcoxon'
        
        
    else
        
        [p1,h1]=ttest(a_s_al,a_ns_al);
        ref_nss=[p1(1) h1(1)];
        
        [p2,h2]=ttest(a_s_al,a_s2);
        ref_180=[p2(1) h2(1)];
        
        [p3,h3]=ttest(a_ns_al,a_ns2);
        ref_ns_180=[p3(1) h3(1)];
        test_a='ttest'
    end
    
    
    res(cr,:)=[ref_nss(2) ref_180(2) ref_ns_180(2)];
    
    metric1=[a_s_al(:,8:12) a_s_al(:,1:7)];
    
    if iii==0
        cl=blushred;
        cl1=squash;
        
    else
        cl=aegean;
        cl1=stone;
    end
    
    f1=figure(1)
    subplot(1,2,cr)
    plot(metric1','Color',cl1)
    ylim([-1 1])
    xlim([1 12])
    xticks([ 3 6 9 12])
    xticklabels({'-90','0','90','180'})
    box('off')
    hold on
    plot(mean(metric1),'k','LineWidth',4,'Color', cl)
    yline(0,'k', 'LineWidth',1,'LineStyle','--')
    ylabel({'Change in tremor severity'})
    xlabel({'Alignment to phase with max';'tremor supression'})
    set(gca,'FontSize',14)
end