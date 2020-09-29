load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/cluster_out_mc.mat','out');
load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','aegean','stone','squash');
color_b1=[stone ; blushred; aegean];
m_ax=1;
low=zeros(10,12);
high=zeros(10,12);
m_arc_a=zeros(10,12);
m_arc_s=zeros(10,12);

cd=1:2; %%%% change here for 1=low amp / 2=high amp
    for iii=1:size(out.start_c,1)
        ref(iii,:)=eval(['nanmedian(out.arc' num2str(cd) '{' num2str(iii) ',2}{m_ax,1})']);
        pre_ns=squeeze(eval(['out.ns_ms{' num2str(iii) ',1}(' num2str(cd) ',:)']));
        ref_ns(1:1e6,1:12)= pre_ns(randi(length(pre_ns),1000000,12));
        [low,high,m_arc_a,m_arc_s]=ns_th(out,ref_ns,iii,low,high,m_arc_a,m_arc_s); clear ref_ns
    end
    % % % %%%% saved output from above
    % % % load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/low_th_group')
    % % % load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/high_th_group')
    
    id_s=[];
    id_a=[];
    for i =1:size(ref,1)
        if (numel(find(ref(i,:)<0)))~=0 % all -except all amplifying subjects
            id_s=[id_s i];
        end
        if(numel(find(ref(i,:)>0)))~=0 % all -except all supressive subjects
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
            
            % main peak at 0 deg
            if phase_peak(ii)==1
                s_al(ii,:)=new_ref(ii,:);
                
            else
                s_al(ii,:)=[new_ref(ii,phase_peak(ii):end) new_ref(ii,1:phase_peak(ii)-1)];
                
                check(ii,1)=sum(diff([phase_peak(ii):size(s_al,2) 1:phase_peak(ii)-1]));
                
            end
            arc{cr,1}=s_al;
        end
        clear s_al
    end
    
    
    low1=low(id_s,:);
    
    % % for iii=1:10
    % % subplot(2,5,iii)
    % % % bar(ref(iii,:))
    % % bar(arc{1,1}(iii,:))
    % % ylim([-1 1])
    % % end
    
    % % figure
    % % subplot(2,1,1)
    % % plot(high','b')
    % % hold on
    % % plot(median(high),'b','lineWidth',4)
    % % plot(arc{1,1}','r')
    % % plot(median(arc{1,1}),'r','lineWidth',4)
    % % subplot(2,1,2)
    % % plot(low1','b')
    % % hold on
    % % plot(median(low1),'b','lineWidth',4)
    % % plot(arc{2,1}','r')
    % % plot(median(arc{2,1}),'r','lineWidth',4)
    
    
    % for g=1:12
    % r(1,g)=kstest(arc{1,1}(:,g)-high(:,g))==1;
    % end
    % sum(r)
    
    for g=1:12
        [j,h]=signrank(arc{1,1}(:,g),(high(:,g)));
        stats.group_amp(g,2)=j;
        stats.group_amp(g,1)=h;clear j h
        %     stats.group_amp(g,1)=h<(0.05/12); clear j h
        [j,h]=signrank(arc{2,1}(:,g),low1(:,g));
        stats.group_sup(g,2)=j;
        stats.group_sup(g,1)=h;clear j h
        %     stats.group_sup(g,1)=h<(0.05/12); clear j h
        dif_a(:,g)=arc{1,1}(:,g)-high(:,g);
        dif_s(:,g)=arc{2,1}(:,g)-low1(:,g);
    end
    
    
    
    load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','aegean','stone','squash','sapphire','azure');
    
    if cd==1
        pp=[sapphire; stone];
    else
        pp=[azure; stone];
    end
    
    f1=figure(50+cd)
    subplot(1,2,1)
    b=bar([median(arc{1,1}(:,1))  median(high(:,1))],'FaceColor','flat','FaceAlpha',.7,'EdgeColor','none','BarWidth',1);
    b.CData(1,:) = pp(1,:);
    b.CData(2,:) = pp(2,:);
    hold on
    plot([arc{1,1}(:,1)' ;  high(:,1)'],'k.')
    ylim([-1 1])
    xticklabels({'stim','non-stim'})
    ylabel({'maximum amplification' ; 'of tremor'})
    box('off')
    set(gca,'FontSize',12);
    
    subplot(1,2,2)
    b=bar([median(arc{2,1}(:,1)) median(low1(:,1))],'FaceColor','flat','FaceAlpha',.7,'EdgeColor','none','BarWidth',1);
    b.CData(1,:) = pp(1,:);
    b.CData(2,:) = pp(2,:);
    hold on
    plot([arc{2,1}(:,1)' ;  low1(:,1)'],'k.')
    ylim([-1 1])
    xticklabels({'stim','non-stim'})
    ylabel({'maximum suppression' ; 'of tremor'})
    box('off')
    set(f1,'color','w');
    set(gca,'FontSize',12);
    