clear all
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\probe SUA_act_mat')
load ('SNR_cycle_probe.mat')
same_a=[28;29];
WaveData_DC=ctx(same_a,:);

hp=find((cellfun(@isempty,SUA))==0);
r=[];
for i=1:length(same_a)
r=[r find(hp(:,1)==same_a(i))];
end

[b,a]=butter(2,[20/(0.5*srn) 30/(0.5*srn)],'bandpass');
            

for i =1:size(ctx(r(1):r(end),:),1)
    psd(i,:)=pwelch(WaveData_DC(i,:),1000,[],1000,1000);
    ecog_filt(i,:)=filtfilt(b,a,WaveData_DC(i,:));
    ecog_env(i,:)=abs(hilbert(ecog_filt(i,:)));
    env=ecog_env(i,:);Ecogfiltered=ecog_filt(i,:);
    onset_raw{i,1}=bursts(env);
    offset_raw{i,1}=bursts_off(env);
    onset_phase_al{i,1}=bursts_aligned(env,Ecogfiltered);
    offset_phase_al{i,1}=bursts_aligned_off(env,Ecogfiltered);
    ecog_phase(i,:)=wrapTo2Pi(angle(hilbert(ecogbf_match(i,:))));
    cy_bursts{i,1}=cycles_10(env,Ecogfiltered);
    clearvars env Ecogfiltered
end


for u=r(1):r(end)
    units_match1=units_match{u,1};
%     close all
    for um=1:size(units_match1,1)
        block = cy_bursts{u-12,1}{2,1}(any(cy_bursts{u-12,1}{2,1},2),:);
        for d1=1:size(block,1)
            for d2=1:size(block,2)
                if d2+1<length(block(d1,:))
                    epoch=block(d1,d2):block(d1,d2+1);
                    l=find(units_match1(um,epoch)==1);
                    pha_b{um,1}{d1,d2}=ecog_phase(u-12,epoch(l));
                    pha_b_all{u-12,um}{d1,d2}=ecog_phase(u-12,epoch(l));
                    units_b_all{u-12,um}{d1,d2}=units_match1(u-12,epoch);
                    idx_spkcycle{um,1}{d1,d2}=epoch(l);
                end
            end
        end
    end
    for cell=1:size(pha_b,1)
        for ii =1:size(pha_b{cell,1},2)
            for i=1:size(pha_b{cell,1},1)
                if ~isempty (pha_b{cell,1}{i,ii})
                    bu{i,1}=pha_b{cell,1}{i,ii}(1);
                end
            end
            bu=bu(~cellfun('isempty',bu));
            vec_lg(cell,ii)=circ_r(cell2mat(bu));
            pref_pha(cell,ii)=circ_mean(cell2mat(bu));
            %         zm(cell,ii) = circ_r(cell2mat(bu)).*(exp(sqrt(-1).*(circ_mean(cell2mat(bu)))));
            clear bu
        end
        cyc_avg(1,:)=mean(vec_lg,1);
        cyc_ang_avg(1,:)=circ_mean(pref_pha);
        
        clearvars -except units_b_all zm cell pha_b_all pha_b vec_lg pl cyc_avg cyc_ang_avg pref_pha units_match cy_bursts ecog_phase u f units_match1
    end
        %per animal
%         for i=11:15
%             if i==11
%                 p1=polarplot([0 cyc_ang_avg(i)], [0, cyc_avg(i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
%                 hold on
%             elseif i==12
%                 p2=polarplot([cyc_ang_avg(i-1) cyc_ang_avg(i)], [cyc_avg(i-1) cyc_avg(i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
%             elseif i==13
%                 p3=polarplot([cyc_ang_avg(i-1) cyc_ang_avg(i)], [cyc_avg(i-1) cyc_avg(i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
%             elseif i==14
%                 p4=polarplot([cyc_ang_avg(i-1) cyc_ang_avg(i)], [cyc_avg(i-1) cyc_avg(i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
%             elseif i==15
%                 p5=polarplot([cyc_ang_avg(i-1) cyc_ang_avg(i)], [cyc_avg(i-1) cyc_avg(i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
%     
%             end
%             rlim([0 0.6])
%         end
%         legend([p1 p2 p3 p4 p5],{'1st cycle','2nd cycle','3rd cycle','4th cycle','5th cycle'},'Box','off','Orientation','horizontal','Location','south','FontSize',9)
%    close all
    clear pha_b
end



for i =1:49
plot(units_b_all{1,1}{i,1})
hold on
end
