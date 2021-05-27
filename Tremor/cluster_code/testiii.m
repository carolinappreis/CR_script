% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/cluster_out_mc.mat','out');
stats=struct;

for cd=1:2 %%%% 1=low amp / 2=high amp
    clearvars -except out cd  stats
m_ax=1;
low=zeros(10,12);
high=zeros(10,12);
m_arc_a=zeros(10,12);
m_arc_s=zeros(10,12);

    for iii=1:size(out.start_c,1)
        ref(iii,:)=eval(['nanmedian(out.arc' num2str(cd) '{' num2str(iii) ',2}{m_ax,1})']); %
        pre_ns=squeeze(eval(['out.ns_ms{' num2str(iii) ',1}(' num2str(cd) ',:)']));
        ref_ns(1:1e6,1:12)= pre_ns(randi(length(pre_ns),1000000,12)); %generate surrogate of non-stim ARC
        [low,high,m_arc_a,m_arc_s]=ns_th(out,ref_ns,iii,low,high,m_arc_a,m_arc_s); clear ref_ns %% align non-stim ARCs to max sup and amp and take 5% and 95% percentile respectivley
    end
    % %% saved output of ns_th for speed
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/low_th_group')
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/high_th_group')

    id_s=[];
    id_a=[];
    for i =1:size(ref,1)
        if (numel(find(ref(i,:)<0)))~=0 % all ARC's -except ARCs that are only amplifying 
            id_s=[id_s i];
        end
        if(numel(find(ref(i,:)>0)))~=0 % all ARC's -except ARCs that are only supressive 
            id_a=[id_a i];
        end
    end
    
    for cr=1:2
        if cr==1
            new_ref=ref(id_a,:);  %  amplifying effect;
            [p,q]=(max(new_ref'));
            phase_peak=q;
        else
            new_ref=ref(id_s,:);  %  supressive effect;
            [p,q]=(min(new_ref'));
            phase_peak=q;
        end
        
        for ii=1:size(new_ref,1)
            
            %%% alignement to max stim effect
            if phase_peak(ii)==1 % main peak at 0 deg
                s_al(ii,:)=new_ref(ii,:);
                
            else
                s_al(ii,:)=[new_ref(ii,phase_peak(ii):end) new_ref(ii,1:phase_peak(ii)-1)];
                
                check(ii,1)=sum(diff([phase_peak(ii):size(s_al,2) 1:phase_peak(ii)-1]));
                
            end
            arc{cr,1}=s_al;
        end
        clear s_al
    end

    %%% use extreme percentiles vs median for stats
    
    mini= high;
    maxi= low;
    
    %median
    mini= m_arc_s;
    maxi= m_arc_a;
    
    
    mini1=mini(id_s,:); %%% match non-stim cases to those where stim ARCs are not exclusivley amplifying --- specific to this data
    
    load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','aegean','stone','squash','sapphire','azure');
    
    if cd==1
        pp=[sapphire; stone];
    else
        pp=[azure; stone];
    end
    
    f1=figure(50+cd)
    %%% plot amplifying effects
    subplot(1,2,1)
    b=bar([median(arc{1,1}(:,1))  median(maxi(:,1))],'FaceColor','flat','FaceAlpha',.7,'EdgeColor','none','BarWidth',1);
    b.CData(1,:) = pp(1,:);
    b.CData(2,:) = pp(2,:);
    hold on
    plot([arc{1,1}(:,1)' ;  maxi(:,1)'],'k.')
    ylim([-1 1])
    xticklabels({'stim','non-stim'})
    ylabel({'maximum amplification' ; 'of tremor'})
    box('off')
    set(gca,'FontSize',12);
    
    %%% plot supressive effects
    subplot(1,2,2)
    b=bar([median(arc{2,1}(:,1)) median(mini1(:,1))],'FaceColor','flat','FaceAlpha',.7,'EdgeColor','none','BarWidth',1);
    b.CData(1,:) = pp(1,:);
    b.CData(2,:) = pp(2,:);
    hold on
    plot([arc{2,1}(:,1)' ;  mini1(:,1)'],'k.')
    ylim([-1 1])
    xticklabels({'stim','non-stim'})
    ylabel({'maximum suppression' ; 'of tremor'})
    box('off')
    set(f1,'color','w');
    set(gca,'FontSize',12);
    f1.Units = 'centimeters';
    f1.OuterPosition= [10, 10, 22, 12];
  
    %%% stats
    for g=1:12
        [j,h]=signrank(arc{1,1}(:,g),(maxi(:,g)));
        stats.group_amp(cd,g,2)=j;
        stats.group_amp(cd,g,1)=h;clear j h
        %     stats.group_amp(g,1)=h<(0.05/12); clear j h
        [j,h]=signrank(arc{2,1}(:,g),mini1(:,g));
        stats.group_sup(cd,g,2)=j;
        stats.group_sup(cd,g,1)=h;clear j h
        %     stats.group_sup(g,1)=h<(0.05/12); clear j h
    end
    
end