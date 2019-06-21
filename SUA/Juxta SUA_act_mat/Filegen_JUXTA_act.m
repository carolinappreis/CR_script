clear all
%cd ('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\Juxta SUA_act_mat\mat')
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/Juxta SUA_act_mat/mat')
load 'animal_region.mat' % comes from code pre_filegen_SUA_act
% load 'non_repeat_animal_region.mat'
region=thal_VA;
for iii=1:size(region,2)
clearvars -except angles_all angles_surr sua_all env_betall ISI_inb ISI_outb spikerate_all idx_spikesin idx_spikesout angles_inb angles_outb nr_bursts region iii A ISI_all spikerate_inbursts spikerate_outbursts ISI_inbursts ISI_outbursts beta_allctx angles_inb idx_spikesur
name=A(region(iii),1:(find(A((region((iii))),:)=='.')-1));
cd ('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\Juxta SUA_act_mat\mat')
load(name);

B=who;

ctxchan=[];
for ii=1:size(B,1)
if ~isempty(min(find(B{ii}=='E')) & min(find(B{ii}=='E')) & min(find(B{ii}=='G')))
ctxchan=ii;
end
end


eval(['samprateold=1/' B{ctxchan} '.interval;']);
eval(['WaveData(1,:)=' B{ctxchan} '.values;']);
WaveData=double(WaveData);
ts=timeseries(WaveData(1,:),0:(1/samprateold):((size(WaveData,2)-1)/samprateold));
ts1=resample(ts,0:0.001:((size(WaveData,2)-1)/samprateold),'linear');
WaveData_DC(1,:)=ts1.data;

cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\Juxta SUA_act_mat\code')
run 'downsamp_bursts.m'
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\Juxta SUA_act_mat\code')
% run 'sprate.m'
run 'sprate_opt.m'

nr_bursts(iii,:)=numel(onset1);
ISI_all{:,iii}=ISI;
spikerate_all{:,iii}=spikerate;
% spikerate_inbursts{:,iii}=[spikes_inb.*(srn/epoch)];
% spikerate_outbursts{:,iii}=[spikes_outb.*(srn/epoch)];
spikerate_inbursts{:,iii}=[spkr_inb];
spikerate_outbursts{:,iii}=[spkr_outb];
ISI_inbursts{:,iii}=ISI_inb;
ISI_outbursts{:,iii}=ISI_outb;
beta_allctx(iii,:)=Ecogfiltered;
env_betall(iii,:)=env;
sua_all{iii,:}=data_ones;
% idx_spikesin{iii,:}=idx_firingin;
% idx_spikesout{iii,:}=idx_firingout;
% idx_spikesur{ii,:}=idx_firingsur;
angles_inb{:,iii}=spkang_inb;
angles_outb{:,iii}=spkang_outb;
angles_surr{:,iii}=spkang_sur;
% nspike_binall{iii,:}=nspike_bin;
angles_all{:,iii}=ang_ones;
end

clearvars -except angles_all angles_surr sua_all tt env_betall  idx_spikesur ISI_inbursts ISI_outbursts spikerate_all nspike_binall idx_spikesin idx_spikesout srn time nr_bursts beta_allctx ISI_all spikerate_inbursts spikerate_outbursts ISI_inbursts ISI_outbursts  angles_inb angles_outb

figure()
cell2mat(ISI_all)
nonzeros(ans)
histogram(ans)
title ('ISI')
% 
figure()
subplot(3,1,1)
histogram(nonzeros(cell2mat(spikerate_all)))
title ('firing rate all time-series')
subplot(3,1,2)
histogram(nonzeros(cell2mat(spikerate_inbursts)))
% boxplot(nonzeros(cell2mat(spikerate_inbursts)))
title ('firing rate in bursts')
subplot(3,1,3)
histogram(nonzeros(cell2mat(spikerate_outbursts)))
title ('firing rate out bursts')

figure ()
subplot(2,1,1)
polarhistogram(cell2mat(angles_inb),12)
title ('phase dist in bursts')
subplot(2,1,2)
%polarhistogram((cell2mat(angles_outb)),12)
polarhistogram(cell2mat(angles_outb),12)
title ('phase dist out bursts')

psi=[abs(mean(exp(1i.*(cell2mat(angles_inb))))) abs(mean(exp(1i.*(cell2mat(angles_surr))))) abs(mean(exp(1i.*(cell2mat(angles_all)))))]

bar(psi)


% if angle(exp(1i.*(m_angprobes)))>0
%     (angle(exp(1i.*(m_angprobes))).*180)./pi
% elseif angle(exp(1i.*(m_angprobes)))<0
%     360-((angle(exp(1i.*(m_angprobes))).*180)./-pi)
% end
% 
%  
% figure()
% plot(time,beta_allctx(1,:))
% hold on
% plot(time(idx_spikesin{1,:}), beta_allctx(1,idx_spikesin{1,:}),'r*')
% plot(time(idx_spikesout{1,:}), beta_allctx(1,idx_spikesout{1,:}),'k.')
% plot(time(sua_all{1,:}), beta_allctx(1,sua_all{1,:}),'bo')
% plot(time,env_betall(1,:))
% plot(time,tt)
% 

