clear all
% cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/Juxta SUA_act_mat/mat')
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\Juxta SUA_act_mat\mat')
load ('SUA_CZ')
load('CZ_cycle')

for i =1:size(stat_d,1)
    UCZ.psd(i,:)=pwelch(WaveData_DC(stat_d(i),:),1000,[],1000,1000);
    UCZ.ecog_filt(i,:)=ecogbf_match(i,:);
    UCZ.ecog_env(i,:)=abs(hilbert(ecogbf_match(i,:)));
    env=UCZ.ecog_env(i,:);Ecogfiltered=UCZ.ecog_filt(i,:);
    UCZ.onset_raw{i,1}=bursts(env);
    UCZ.offset_raw{i,1}=bursts_off(env);
    UCZ.onset_phase_al{i,1}=bursts_aligned(env,Ecogfiltered);
    UCZ.offset_phase_al{i,1}=bursts_aligned_off(env,Ecogfiltered);
    UCZ.ecog_phase(i,:)=wrapTo2Pi(angle(hilbert(ecogbf_match(i,:))));
    cy_bursts{i,1}=cycles_10(env,Ecogfiltered);
    clearvars env Ecogfiltered
end

for f=1:size(cy_bursts,1)
    block = cy_bursts{f,1}{2,1}(any(cy_bursts{f,1}{2,1},2),:);
    for d1=1:size(block,1)
        for d2=1:size(block,2)
            if d2+1<length(block(d1,:))
                epoch=block(d1,d2):block(d1,d2+1);
                l=find(units_match(f,epoch)==1);
                pha_b{f,1}{d1,d2}=UCZ.ecog_phase(f,epoch(l));
                idx_spkcycle{f,1}{d1,d2}=epoch(l);
            end
        end
    end
end


for cell=1:size(pha_b,1)
    for ii =1:size(pha_b{cell,1},2)
        for i=1:size(pha_b{cell,1},1)
            if ~isempty (pha_b{cell,1}{i,ii})
                    bu{i,1}=mean(exp(sqrt(-1).*(pha_b{cell,1}{i,ii})));        
            end
        end
        bu=bu(~cellfun('isempty',bu));
        vec_lg(cell,ii)=abs(mean(cell2mat(bu)));
        pref_pha(cell,ii)=angle(mean(cell2mat(bu)));
        zm(cell,ii) = (abs(mean(cell2mat(bu)))).* (mean(cell2mat(bu)));
        clear bu
    end
    cyc_avg(1,:)=mean(vec_lg,1);
    cyc_ang_avg(1,:)=circ_mean(pref_pha);
    
    clearvars -except zm cell pha_b vec_lg  cyc_avg cyc_ang_avg pref_pha
end



for i=11:15
    
  
    if i==11
        p1=polarplot([0 cyc_ang_avg(i)], [0, cyc_avg(i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
          hold on
    elseif i==12
        p2=polarplot([cyc_ang_avg(i-1) cyc_ang_avg(i)], [cyc_avg(i-1) cyc_avg(i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
    elseif i==13
        p3=polarplot([cyc_ang_avg(i-1) cyc_ang_avg(i)], [cyc_avg(i-1) cyc_avg(i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
    elseif i==14
        p4=polarplot([cyc_ang_avg(i-1) cyc_ang_avg(i)], [cyc_avg(i-1) cyc_avg(i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
    elseif i==15
        p5=polarplot([cyc_ang_avg(i-1) cyc_ang_avg(i)], [cyc_avg(i-1) cyc_avg(i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
        
    end
    rlim([0 0.6])
end
    legend([p1 p2 p3 p4 p5],{'1st cycle','2nd cycle','3rd cycle','4th cycle','5th cycle'},'Box','off','Orientation','horizontal','Location','south','FontSize',9)



for n=1:size(pha_b,1)
    for i=11:15
%         subplot(1,size(pha_b,1),n)
figure(n+1)
        if i==11
            p1=polarplot([0 pref_pha(n,i)], [0, vec_lg(n,i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
                    hold on
        elseif i==12
            p2=polarplot([pref_pha(n,i-1) pref_pha(n,i)], [vec_lg(n,i-1), vec_lg(n,i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
        elseif i==13
            p3=polarplot([pref_pha(n,i-1) pref_pha(n,i)], [vec_lg(n,i-1), vec_lg(n,i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
        elseif i==14
            p4=polarplot([pref_pha(n,i-1) pref_pha(n,i)], [vec_lg(n,i-1), vec_lg(n,i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
        elseif i==15
            p5=polarplot([pref_pha(n,i-1) pref_pha(n,i)], [vec_lg(n,i-1), vec_lg(n,i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
            
        end
    end
    rlim([0 0.7])
    legend([p1 p2 p3 p4 p5],{'1st cycle','2nd cycle','3rd cycle','4th cycle','5th cycle'},'Box','off','Orientation','horizontal','Location','south','FontSize',9)
end

