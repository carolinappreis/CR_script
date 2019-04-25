clear all
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
 cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A3_Thal/mat')
load ('data_all' , 'freq')

 cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/Juxta SUA_act_mat/mat')
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\Juxta SUA_act_mat\mat')
load 'animal_region.mat' % comes from code pre_filegen_SUA_act
% load 'non_repeat_animal_region.mat'
% load 'single_subj_VM_VL.mat'
% all_regions={thal_VL; thal_VM};
thal_CZ= thal_VL;
thal_BZ=[thal_VM thal_VA];
% all_regions={thal_CZ; thal_BZ;}
all_regions={thal_BZ;}
freq=[{1:7} {8:14} {15:35} {36:80} {81:100} ];

for t=1:size(freq,2)
    for i=1:size(all_regions,1)
        for  j=1:size(all_regions{i,:},2)
            
            clearvars -except A all_regions i j ang ang_cons regional_ang m data_all euler euler1 ang_mean ang_cons t freq ang1
            name=A(all_regions{i,:}(j),1:(find(A((all_regions{i,:}(j)),:)=='.')-1));
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
            
            sr=1/unite.interval;
            srn=1000;
            dataold=unite.values';
            dataold=full(dataold);
            data=zeros(1,100000);
            timeold=0:1/sr:(size(dataold,2)-1)/sr;
            time=0:1/srn:(size(data,2)-1)/srn;
            spk_t=timeold(find(dataold==1));
            spk_tround=round(spk_t,3);
            nn=[];
            for ii=1:length(spk_t)
                [ d, ix ] = min( abs( time-spk_tround(ii) ) );
                nn=[nn ix];
            end
            data(nn)=1;
            data_all{i,j}=data;
            data_ones=find(data==1);
             [b,a]=butter(2,[freq{1,t}(1)/(0.5*srn) freq{1,t}(end)/(0.5*srn)],'bandpass');
            Ecogfiltered=filtfilt(b,a,WaveData_DC);
            ang{t,i,:}{:,j}=angle(hilbert(Ecogfiltered(data_ones)));
            euler{t,i,:}(1,j)=(sum(exp(sqrt(-1)*(ang{t,i,:}{:,j})))./(length(ang{t,i,:}{:,j})));
        end
        euler1(t,i)=mean(euler{t,i,:});
        ang1{t,i}=cell2mat(ang{t,i});
    end
end

bar(abs(euler1),'FaceColor',[0.5 0 0.5])
title ('BZ')
 xticklabels({'1-7','8-12','15-35','36-80','81-100'})
ylabel ('Spike-EcoG-{\beta} coherence')
xlabel ('Frequency (Hz)')
box('off')

% cd('/Users/Carolina/Desktop/TC_data') 
% savefig('sua_phaseconsist')
% saveas(gca,'sua_phaseconsist.png')

% for t=1:4
%     subplot (1,4,t)
%     polarhistogram(ang1{t,3},12)
% end
% title ('mean phase VM')

% figure()
% for i=1:size(data_all,2)
%     subplot(size(data_all,2),1,i)
%     plot(time,data_all{3,i})
%     ylim ([-2 2])
% %     xlim ([0 10])
% end


% plot (time,Ecogfiltered)
% xlim ([0 0.1])
% figure()
% plot(time,data_all{1,1})
% xlim ([0 0.1])