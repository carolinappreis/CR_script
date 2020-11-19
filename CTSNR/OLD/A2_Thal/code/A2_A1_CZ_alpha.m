clear all
cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A3_Thal/mat')
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
load CZ_ctx_probe.mat
freq=[5:15];
samprate=1000;
n=1;
for subj=1:size(data,1)
    for r=2:size(data{subj,1},1)
        [Pxx_indA,F_ind]=mscohere(data{subj,1}(r,:),data{subj,1}(1,:),samprate,[],samprate,samprate); %Magnitude-squared coherence between ctx-a given contact
        if sum(Pxx_indA(freq))./sum(Pxx_indA(1:end))>0.1
            p(n,:)=sum(Pxx_indA(freq))./sum(Pxx_indA(1:end));
            freq_peak= find(Pxx_indA==(max(Pxx_indA(freq,1))));
            p1(n,:)=freq_peak;
            p2(n,:)=max(Pxx_indA(freq,1));
            [b,a]=butter(2,[(freq_peak-2)/(0.5*samprate) (freq_peak+2)/(0.5*samprate)],'bandpass');
                        [powerA,f]=pwelch(data{subj,1}(r,:),1000,[],1000,1000);
                        plot(log(powerA))
                        xlim([0 50])
                        hold on
            [bb,aa]=butter(2,[15/(0.5*samprate) 35/(0.5*samprate)],'bandpass');
            Ecogfiltered= filtfilt(bb,aa,data{subj,1}(1,:));
            env=abs(hilbert(Ecogfiltered));
            raw_onset=bursts(env);
            raw_onset=double(raw_onset{2,1});
            aligned_onset=bursts_aligned(env,Ecogfiltered);
            aligned_onset=double(aligned_onset{2,1});
            CZ_filt=filtfilt(b,a,data{subj,1}(r,:));
            CZ_env=abs(hilbert(CZ_filt));
            
            el=500;
            for ii=1:length(raw_onset)
                if raw_onset(ii)>el
                    epochs_raw1(ii,:)= (CZ_env(raw_onset(ii)-el:raw_onset(ii)+el)-median(CZ_env(raw_onset(ii)-el:raw_onset(ii))))./(median(CZ_env(raw_onset(ii)-el:raw_onset(ii))));
                    epochs_aligned1(ii,:)= data{subj,1}(r,(aligned_onset(ii)-el:aligned_onset(ii)+el));
                end
            end
            
            epochs_raw(n,:)=median(epochs_raw1);
            epochs_aligned(n,:)=median(epochs_aligned1);
            clear epochs_aligned1 epochs_raw1
            n=n+1;
        end
    end
end

plot(median(epochs_raw))
figure()
plot(sum(epochs_aligned))
xlim([300 700])

