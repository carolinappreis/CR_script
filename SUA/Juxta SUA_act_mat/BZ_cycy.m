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
                bu{i,1}=pha_b{cell,1}{i,ii}(1);      
            end
        end
        bu=bu(~cellfun('isempty',bu));
        vec_lg(cell,ii)=circ_r(cell2mat(bu));
        pref_pha(cell,ii)=circ_mean(cell2mat(bu));
        zm(cell,ii) = circ_r(cell2mat(bu)).*(exp(sqrt(-1).*(circ_mean(cell2mat(bu)))));
        clear bu
    end
    cyc_avg(1,:)=mean(vec_lg,1);
    cyc_ang_avg(1,:)=circ_mean(pref_pha);
    
    clearvars -except zm cell pha_b vec_lg  cyc_avg cyc_ang_avg pref_pha
end


for i=10:14
    zm1=mean(zm,1);
  
    if i==10
        p1=polarplot([0 real(zm1(i))], [0, imag(zm1(i))],'linewidth',2,'MarkerIndices',[2],'Marker','d')
          hold on
    elseif i==11
        p2=polarplot([real(zm1(i-1)) real(zm1(i))], [imag(zm1(i-1)), imag(zm1(i))],'linewidth',2,'MarkerIndices',[2],'Marker','d')
    elseif i==12
        p3=polarplot([real(zm1(i-1)) real(zm1(i))], [imag(zm1(i-1)), imag(zm1(i))],'linewidth',2,'MarkerIndices',[2],'Marker','d')
    elseif i==13
        p4=polarplot([real(zm1(i-1)) real(zm1(i))], [imag(zm1(i-1)), imag(zm1(i))],'linewidth',2,'MarkerIndices',[2],'Marker','d')
    elseif i==14
        p5=polarplot([real(zm1(i-1)) real(zm1(i))], [imag(zm1(i-1)), imag(zm1(i))],'linewidth',2,'MarkerIndices',[2],'Marker','d')
        
    end
    rlim([0 0.6])
end
legend([p1 p2 p3 p4 p5],{'1st cycle','2nd cycle','3rd cycle','4th cycle','5th cycle'})

figure()
for n=1:size(pha_b,1)
    for i=10:14
        subplot(1,size(pha_b,1),n)
        if i==10
            p1=polarplot([0 real(zm(n,i))], [0, imag(zm(n,i))],'linewidth',2,'MarkerIndices',[2],'Marker','d')
                    hold on
        elseif i==11
            p2=polarplot([real(zm(n,i-1)) real(zm(n,i))], [imag(zm(n,i-1)), imag(zm(n,i))],'linewidth',2,'MarkerIndices',[2],'Marker','d')
        elseif i==12
            p3=polarplot([real(zm(n,i-1)) real(zm(n,i))], [imag(zm(n,i-1)), imag(zm(n,i))],'linewidth',2,'MarkerIndices',[2],'Marker','d')
        elseif i==13
            p4=polarplot([real(zm(n,i-1)) real(zm(n,i))], [imag(zm(n,i-1)), imag(zm(n,i))],'linewidth',2,'MarkerIndices',[2],'Marker','d')
        elseif i==14
            p5=polarplot([real(zm(n,i-1)) real(zm(n,i))], [imag(zm(n,i-1)), imag(zm(n,i))],'linewidth',2,'MarkerIndices',[2],'Marker','d')
            
        end
    end
    rlim([0 0.6])
end
legend([p1 p2 p3 p4 p5],{'1st cycle','2nd cycle','3rd cycle','4th cycle','5th cycle'})

% circ_plot(cell2mat(bu),'pretty','bo',true,'linewidth',2,'color','r')

%checking plots
% plot(log(UCZ.psd(:,1:50))')
% env=UCZ.ecog_env(1,:);Ecogfiltered=UCZ.ecog_filt(1,:);

% plot(time,env)
% hold on
% plot(time,Ecogfiltered)
% plot(time(UCZ.onset_phase_al{1,1}{1,1}),Ecogfiltered(UCZ.onset_phase_al{1,1}{1,1}),'bo')
% plot(time(UCZ.offset_phase_al{1,1}{1,1}),Ecogfiltered(UCZ.offset_phase_al{1,1}{1,1}),'ko')
% plot(time(UCZ.onset_raw{1,1}{1,1}),env(UCZ.onset_raw{1,1}{1,1}),'b.')
% plot(time(UCZ.offset_raw{1,1}{1,1}),env(UCZ.offset_raw{1,1}{1,1}),'k.')
%
% plot(time,Ecogfiltered)
% hold on
% plot(time(cy_bursts{1,1}{2,1}(2,:)),Ecogfiltered(cy_bursts{1,1}{2,1}(2,:)),'r.')
% plot(time(cell2mat(idx_spkcycle{1,1}(1,:))),Ecogfiltered(cell2mat(idx_spkcycle{1,1}(1,:))),'k.')
