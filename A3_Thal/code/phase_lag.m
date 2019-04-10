
clear all

Fs=1000;
%cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A2_Thal/mat')
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
load SNR_ctx_probesraw.mat
WaveSNR=WaveData_DCall;
idx_SNR=[];
for i =1:29
    if size(WaveSNR{i,1},1)~=1
        idx_SNR=[idx_SNR i];
    end
end



%cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A2_Thal/mat')
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
load BZ_ctx_probesraw.mat
WaveBZ=WaveData_DCall;
idx_BZ=[];
for i =1:29
    if size(WaveBZ{i,1},1)~=1
        idx_BZ=[idx_BZ i];
    end
end


subj=[idx_BZ(idx_BZ==idx_SNR)];
for i=1:length(subj)
    m=1;
for n=2:(min([ size(WaveBZ{subj(i),1},1)  size(WaveSNR{subj(i),1},1)]));
%         [Cxy,f] = mscohere((WaveBZ{subj(i),1}(n,:)),(WaveSNR{subj(i),1}(n,:)),Fs,[],Fs,Fs);
%         [valP,idxP] = (findpeaks(Cxy(1:100)));
%         poi=sort(valP,'descend');
%         if idxP(find(valP==(poi(1))))>5
%             filtrange=idxP(find(valP==(poi(1))));
%         else
%             filtrange=idxP(find(valP==(poi(2))));
%         end
        filtrange=25;
        [b,a]=butter(2,[(filtrange-5)/(0.5*Fs) (filtrange+5)/(0.5*Fs)],'bandpass');
        filtsign(1,:)=filtfilt(b,a,WaveBZ{subj(i),1}(n,:));
        filtsign(2,:)=filtfilt(b,a,WaveSNR{subj(i),1}(n,:));
        filtsign(3,:)=filtfilt(b,a,WaveData_DCall{subj(i),1}(1,:));
        dif_angs(i,m,:)=angle(hilbert(filtsign(1,:)))-angle(hilbert(filtsign(2,:)));
        m=m+1;
        
    end
    
end

for i=1:size(dif_angs,1)
    for j=1:size(dif_angs,2)
        if dif_angs(i,j)>pi
            dif_angs(i,j)=dif_angs(i,j)-2*pi;
        elseif dif_angs(i,j)<(-pi)
            dif_angs(i,j)=dif_angs(i,j)+2*pi;
        end
    end
end

for i=1:size(dif_angs,1)
psi(i,:)=abs(mean(exp(1i*(dif_angs(i,:)))))
end


polar(repmat((dif_angs(6,:)),1,2)',repmat([0 1],1, 100000)');
hold on
polar(repmat(abs(mean(exp(1i*(dif_angs(6,:))))),1,2)',repmat([0 1],1, 1)');

polarhistogram(dif_angs(6,:),12);

