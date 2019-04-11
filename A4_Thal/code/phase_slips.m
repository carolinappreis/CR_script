clear all
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal')

load('SNR.mat');
% for ii=1:length(SNR.idrat)
% data{ii,1}=SNR.filt_thal{SNR.idrat(ii),1}
% end
% data=vertcat(data{:});
SNR.idrat=[1:14];
for ik=1:length(SNR.idrat)
    for ct=1:size(SNR.phase_thal{SNR.idrat(ik),1},1)
        non_nomr=[];epochs_z=[];epochs_z1=[];
        ref1=SNR.onset_raw{1,(SNR.idrat(ik))}{2,1};
        ref2=SNR.offset_raw{1,SNR.idrat(ik)}{2,1};
        [dur,dur_idx]=sort(ref2-ref1,'ascend');
        ref3=ref2(dur_idx(end-24:end));
        non_norm=SNR.phase_ctx(SNR.idrat(ik),:)-SNR.phase_thal{SNR.idrat(ik),1}(ct,:);
        for x =1:size(non_norm,2)
            if non_norm(1,x)>pi
                non_norm(1,x)=(non_norm(1,x))-(2.*pi);
            elseif non_norm(1,x)<-pi
                non_norm(1,x)=(non_norm(1,x))+(2.*pi);
            else
                non_norm(1,x)= non_norm(1,x);
            end
        end
        
        znon_norm=zscore(non_norm);
        el=400;
        for ii=1:length(ref3)
            if ref3(ii)>el
                epochs_z(ii,:)=diff(znon_norm(ref3(ii)-el:ref3(ii)+el));
            end
            for ii=1:size(epochs_z,1)
                for ff=1:length(epochs_z(ii,:))
                    if epochs_z(ii,ff)>1,96 | epochs_z(ii,ff)<-1,96;
                        epochs_z1(ii,ff)=1;
                    else
                        epochs_z1(ii,ff)=0;
                    end
                end
            end
        end
        %         plot(epochs_z1(ii,:),'r.')
        %         imagesc(epochs_z1)
        %         xticks([200:200:800])
        %         xlim([200 800])
        %         xticklabels ({'-200','0','200','400'})
        epochs_zd(ik,ct,:,:)=epochs_z1;
    end
end

epochs_zd1=reshape(epochs_zd,size(epochs_zd,1)*size(epochs_zd,2),size(epochs_zd,3),size(epochs_zd,4));
slip_b=squeeze(mean(epochs_zd1,1));

imagesc(slip_b)
xlabel ('Time(msec)')
ylabel('Bursts # (Sorted by length)')
xticks([200:200:800])
xlim([200 800])
xticklabels ({'-200','0','200','400'})
title('SNR')

% if squeeze(sum(epochs_zd1(1:end,1,:))./size(epochs_zd1,3))==slip_b(1,:)'
%     cr=1
% end

% slip_b(slip_b~=0)
