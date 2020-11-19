clear all
close all
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A3_Thal/mat')
load ('data_all' , 'freq')

% cd ('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\Juxta SUA_act_mat\mat')
cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/Juxta SUA_act_mat/mat')
load ('animal_region.mat') % comes from code pre_filegen_SUA_act
regions={thal_VA ;thal_VL;thal_VM};

for xx=1:size(regions,1)
    for t=1:size(freq,1)
        for iii=1:size(regions{xx,1},2)
            % clearvars -except angles_all angles_surr sua_all env_betall ISI_inb ISI_outb spikerate_all idx_spikesin idx_spikesout angles_inb angles_outb nr_bursts region iii A ISI_all spikerate_inbursts spikerate_outbursts ISI_inbursts ISI_outbursts beta_allctx angles_inb idx_spikesur
            name=A(regions{xx,1}(iii),1:(find(A((regions{xx,1}((iii))),:)=='.')-1));
%            cd ('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\Juxta SUA_act_mat\mat')
              cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/Juxta SUA_act_mat/mat')
            load(name);
            B=who;
            ctxchan=[];
            for ii=1:size(B,1)
                if ~isempty(min(find(B{ii}=='E')) & min(find(B{ii}=='E')) & min(find(B{ii}=='G')))
                    ctxchan=ii;
                end
            end
            
            WaveData=[];
            eval(['samprateold=1/' B{ctxchan} '.interval;']);
            eval(['WaveData(1,:)=' B{ctxchan} '.values;']);
            WaveData=double(WaveData);
            ts=timeseries(WaveData(1,:),0:(1/samprateold):((size(WaveData,2)-1)/samprateold));
            ts1=resample(ts,0:0.001:((size(WaveData,2)-1)/samprateold),'linear');
            WaveData_DC(1,:)=ts1.data;
            sr=1/unite.interval;
            srn=1000;
            dataold=unite.values';
            dataold=full(dataold);
            data=zeros(1,100000);
            timeold=0:1/sr:(size(dataold,2)-1)/sr;
            time=0:1/srn:(size(data,2)-1)/srn;
            spk_t=timeold(find(dataold==1));
            spk_tround=round(spk_t,3);
            n=[];
            for i=1:length(spk_t)
                [ d, ix ] = min( abs( time-spk_tround(i) ) );
                n=[n ix];
            end
            data(n)=1;
            data_ones=find(data==1);
            
            [b,a]=butter(2,[(freq(t)-5)/(0.5*srn) (freq(t)+5)/(0.5*srn)],'bandpass');
            if t==length(freq)
                [b,a]=butter(2,[49/(0.5*srn) 60/(0.5*srn)],'bandpass');
            end
            Ecogfiltered=filtfilt(b,a,WaveData_DC);
            env=abs(hilbert(Ecogfiltered));
            ang=angle(hilbert(Ecogfiltered));
            threshold=prctile(env,75);
            tt(size(env,1):size(env,2))=threshold;
            indexexceed=find(env>threshold);
            diffindex=diff(indexexceed);
            pnts=find(diffindex>1);
            begin=indexexceed(pnts+1);
            ending=indexexceed(pnts);
            begin2=[indexexceed(1) begin];
            ending2=[ending indexexceed(end)];
            
            duration=ending2-begin2;
            mean_b=mean(duration);
            
            ind_b=[];
            for i=1:(length(begin2))
                if (ending2(i)-begin2(i))>=100 %original
                    %                 if (ending2(i)-begin2(i))>=mean_b
                    %                 if (ending2(i)-begin2(i))>=50 && (ending2(i)-begin2(i))<=mean_b
                    ind_b=[ind_b i];
                end
            end
            
            begin3=begin2(ind_b);
            ending3=ending2(ind_b);
            
            thr=200;
            ind_b1=[];
            for i=1:(length(begin3)-2)
                if (begin3(i+1)-ending3(i))>=thr && (begin3(i+2)-ending3(i+1))>=thr
                    ind_b1=[ind_b1 i+1];
                end
            end
            
                
            if (begin3(2)-ending(1))>=thr
                ind_b1= [1 ind_b1];
            end
            
            if (begin3(length(begin3))-ending3(length(begin3)-1))>=thr
                ind_b1=[ind_b1 length(begin3)];
            end
            
            onset1=begin3(ind_b1);
            offset1=ending3(ind_b1);
            
            m=1;
            for i =1:srn:(length(data)-srn);
                spikerate_all{iii,xx}(1,m)=numel(find(data(i:i+srn)==1));
                m=m+1;
            end
            
            m=1;
            for j=1:length(onset1)
                output_in{:,j}= (onset1(j)+find(data(onset1(j):offset1(j))==1))-1;
                ISI_inb1{:,j}=diff(time(find(data(onset1(j):offset1(j))==1)));
                spkr_inb(1,m)=numel(find(data(onset1(j):offset1(j))==1)).*srn./(length(time(onset1(j):offset1(j))));
                m=m+1;
            end
            
            for a=1:(length(data_ones))
                idx_sur=randi([thr+1,(length(data)-thr)],1,1);
                output_sur{:,a}=idx_sur+find(data(idx_sur:idx_sur+thr)==1);
                ISI_sur1{:,a}=diff(time(find(data(idx_sur:idx_sur+thr)==1)));
                spkr_sur(1,a)=numel(find(data(idx_sur:idx_sur+thr)==1)).*srn./(length(time(idx_sur:idx_sur+thr)));
            end
            
            spikerate_inbursts{t,iii,xx}=nonzeros(spkr_inb);
            spikerate_sur{t,iii,xx}=nonzeros(spkr_sur);
            
            ISI_all{iii,xx}=diff(time(data==1));
            ISI_inb{t,iii,xx}= nonzeros(cell2mat(ISI_inb1));
            ISI_sur{t,iii,xx}= nonzeros(cell2mat(ISI_sur1));
            
            idx_firingin=nonzeros(cell2mat(output_in));
            idx_firingsur=nonzeros(cell2mat(output_sur));
            spkang_inb{t,iii,xx}=ang(idx_firingin); % ang= phases of filtered beta
            spkang_sur{t,iii,xx}=ang(idx_firingsur);
            spkang_all{t,iii,xx}=ang(data_ones);
            
            clearvars Ecogfiltered output_sur ISI_sur1  spkr_sur output_in ISI_inb1 spkr_inb
        end
    end
end

clearvars -except spkang_all spkang_inb spkang_sur ISI_all ISI_inb ISI_sur spikerate_inbursts spikerate_all spikerate_sur data_all


%%%------------------------------------------------------------------plots

spikerate=[];
for i=1:size(spikerate_all,2)
    for j=1:size(spikerate_all,1)
        spikerate=[spikerate double(spikerate_all{j,i})];
    end
    subplot(3,1,i)
    histogram(spikerate,25)
    xlim ([0 50])
    xticks([0:10:50])
    xlabel ('frequency (Hz)')
    ylabel ('count')
end
% 
ISI=[];
for i=1:size(spikerate_all,2)
    for j=1:size(spikerate_all,1)
        ISI=[ISI double(ISI_all{j,i})];
    end
    subplot(3,1,i)
    histogram(ISI,100)
    %     xlim ([0 50])
    %     xticks([0:10:50])
    xlabel ('Time (sec)')
    ylabel ('count')
    xlim([0 2.5])
end
% 
% ISI per contact
% for j=1:3
%     subplot(3,1,j)
%     histogram(ISI_all{j,1})
%     %     xlim ([0 0.5])
%     xlabel ('Time(sec)')
%     ylabel ('count')
% end
% 
for f=1:size(spikerate_inbursts,3)
    for i=1:size(spikerate_inbursts,1)
        for j=1:size(spikerate_inbursts,2)
            if ~isempty (spikerate_inbursts{i,j,f})
            m_sr_b(i,j)=mean(double(spikerate_inbursts{i,j,f}));
            m_sr_s(i,j)=mean(double(spikerate_sur{i,j,f}));
            end
        end
        stats_sr1(i,:)=ranksum(abs(m_sr_b(i,:)), abs(m_sr_s(i,:)))<0.05/(size(m_sr_b(i,:),2));
    end
    stats_sr(:,f)=stats_sr1;
    skr_inb_region(:,f)=nanmean(m_sr_b,2);
    skr_sur_region(:,f)=nanmean(m_sr_s,2);
end
n=[];
m=1;
for i=1:3
    n(:,m:m+1)=[(skr_inb_region(:,i)) (skr_sur_region(:,i))];
    m=m+2;
end
figure()
bar(n)
xlabel ('Cortical signal')
ylabel ('Firing rate(Hz)')
xticklabels({'5-15Hz','16-26Hz','27-37Hz','38-48Hz','49-100Hz'})
title('Firing rate in bursts vs. surr')

% 

for f=1:3
    for i=1:5
        for j=1:7
            if ~isempty (spkang_inb{i,j,f})
                m_a_b(i,j)=mean(exp(sqrt(-1)*(double(spkang_inb{i,j,f}))));
                m_a_s(i,j)=mean(exp(sqrt(-1)*(double(spkang_sur{i,j,f}))));
            end
        end
        stats(i,f)=ranksum(abs(m_a_b(i,:)), abs(m_a_s(i,:)));
        stats_psi(i,f)=ranksum(abs(m_a_b(i,:)), abs(m_a_s(i,:)))<0.05/(size(m_a_b(i,:),2));
    end
    ska_inb_region(:,f)=abs(mean(m_a_b,2));
    ska_sur_region(:,f)=abs(mean(m_a_s,2));
end

n=[];
m=1;
for i=1:3
    n(:,m:m+1)=[(ska_inb_region(:,i)) (ska_sur_region(:,i))];
    m=m+2;
end
figure()
bar(n)
xlabel ('Cortical signal')
ylabel ('PSI')
xticklabels({'5-15Hz','16-26Hz','27-37Hz','38-48Hz','49-100Hz'})
title('PSI in bursts vs surrogates')

% cd('/Users/Carolina/Desktop/TC_data') 
% savefig('sua_frate')
% saveas(gca,'sua_frate.png')
% 
% cd('/Users/Carolina/Desktop/TC_data') 
% savefig('sua_phaseconsist')
% saveas(gca,'sua_phaseconsist.png')

% bar(rad2deg(angle(n)))
% xlabel ('Cortical signal')
% ylabel ('Angle (deg)')
% xticklabels({'5-15Hz','16-26Hz','27-37Hz','38-48Hz'})
% title('Mean angle INSIDE BURSTS Juxta')
%
%
% m_angprobes=(rad2deg(angle(ska_inb_region)));
% for i=1:5
%     for j=1:3
%         if angle(exp(1i.*(m_angprobes(i,j))))>0
%             m_angprobes(i,j)=(angle(exp(1i.*(m_angprobes(i,j)))).*180)./pi
%         elseif angle(exp(1i.*(m_angprobes(i,j))))<0
%             m_angprobes(i,j)=360-((angle(exp(1i.*(m_angprobes(i,j)))).*180)./-pi)
%         end
%     end
% end

% bar(m_angprobes)



