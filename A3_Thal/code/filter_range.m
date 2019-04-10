
clear all
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A3_Thal/mat')
load 'data_all.mat'

samprate=1000;
m=1;
freq_ms={5:15; 16:26 ; 27:37 ;38:48; 49:100};
for l=1:size(freq,1);
    for i=1:size(data_all,1)
        for r=2:size(data_all{i,:},1)
            %         [power,f]=pwelch(data_all(r,:),1000,[],1000,1000);
            [Pxx_ind,F_ind]=mscohere(data_all{i,:}(1,:),data_all{i,:}(r,:),samprate,[],samprate,samprate);
            Pxx_ind1=Pxx_ind(freq_ms{l});
            filtrange{l,1}(m,:)=(freq_ms{l}(1,1))+find(Pxx_ind1==max(Pxx_ind1));
            m=m+1;
        end
    end
    freq(l,:)=round(mean(nonzeros(filtrange{l,1})))
end

clearvars -except data_all freq
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A3_Thal/mat')
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
freq=([alpha_mean betal_mean betah_mean gamma_mean]);
% save 'data_all'
clear all
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A3_Thal/mat')
load ('data_all.mat')

xticklabels({'5-15Hz','16-26Hz','27-37Hz','38-48Hz','49-100Hz'})