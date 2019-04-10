
clear all
close all
Fs=1000;
cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A3_Thal/mat')
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
load BZ_ctx_probe.mat
WaveA=WaveData_DCall;
dataA=data;
idx_A=[];
for i =1:29
    if size(WaveA{i,1},1)~=1
        idx_A=[idx_A i];
    end
end

cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A3_Thal/mat')
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
load CZ_ctx_probe.mat
WaveB=WaveData_DCall;
dataB=data;
idx_B=[];
for i =1:29
    if size(WaveB{i,1},1)~=1
        idx_B=[idx_B i];
    end
end

load ('data_all.mat','freq')
subj= idx_B(ismember(idx_B,idx_A));
samprate=1000;

for t=1:size(freq,1)
    psi_rec=[];
    ctcts=0;
    for r=1:length(subj)
        coh_thal=[];
        clearvars a b euler
        pp=1;
        
        if t==length(freq)
            [b,a]=butter(2,[49/(0.5*samprate) 100/(0.5*samprate)],'bandpass');
        else
            [b,a]=butter(2,[(freq(t)-5)/(0.5*samprate) (freq(t)+5)/(0.5*samprate)],'bandpass');
        end
        for ra=2:size(WaveA{subj(r),1},1);
            for rb=2:size(WaveB{subj(r),1},1);
                
                [Pxx_ind,F_ind]=mscohere(WaveB{subj(r),1}(rb,:),WaveA{subj(r),1}(ra,:),samprate,[],samprate,samprate); %Magnitude-squared coherence between ctx-a given contact            
                if  sum(Pxx_ind(freq(t)-5:freq(t)+5))./sum(Pxx_ind(1:end))>0.1
                    ctcts=ctcts+1;
                    coh_thal=[coh_thal sum(Pxx_ind(freq(t)-5:freq(t)+5))./sum(Pxx_ind(1:end))];
                    
                    non_norm=angle(hilbert(filtfilt(b,a,WaveB{subj(r),1}(rb,:))))-angle(hilbert(filtfilt(b,a,WaveA{subj(r),1}(ra,:))));
                    
                    for x =1:size(non_norm,2)
                        if non_norm(1,x)>pi
                            non_norm(1,x)=(non_norm(1,x))-(2.*pi);
                        elseif non_norm(1,x)<-pi
                            non_norm(1,x)=(non_norm(1,x))+(2.*pi);
                        else
                            non_norm(1,x)= non_norm(1,x);
                        end
                    end
                    
                    dif_angs=non_norm;
                    clearvars non_norm;
                    euler1=abs(sum(exp(sqrt(-1)*(dif_angs)))./(length(dif_angs)));
                    euler(pp,:)=euler1;
                    pp=pp+1;
                    
                    psi_rec(r,:)=mean(euler);
                    clearvars dif_angs
                    
                end
            end
        end
    end
    if ~isempty (psi_rec)
        psi_all(t)=mean(psi_rec);
    else
        psi_all(t)=NaN;
    end
    ctcts_all(t,1)=ctcts;

end
bar(psi_all,0.5)
xticklabels({'5-15','16-26','27-37','38-48','49-100'})
xlabel ('Frequencies(Hz)')
ylabel ('PSI BZ-SNR')
box off

cd('/Users/Carolina/Desktop/TC_data') 
savefig('BUA_BZ_RT')
saveas(gca,'BUA_BZ_RT.png')

% 
% 
% 
% figure()
% for i=1:length(freq);
%     title('PSI')
%     subplot(length(freq),1,i)
%     imagesc(psi{i,1});
%     colorbar
% end
% %
% for i=1:2;
%     figure(2)
%     subplot(2,1,i)
%     imagesc(squeeze(ang(i,:,:))');
% end
% %
% test=squeeze(dif_angs(2,4,2,:));
% figure()
% subplot(2,1,1)
% polar(repmat(test,2,1)',repmat([0 1],1, 100000));
% subplot(2,1,2)
% polarhistogram(test)
