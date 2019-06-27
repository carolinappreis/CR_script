clear all
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal')
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal')
load('BZ_opt.mat');

for duration=1:2
for ik=1:size(BZ.filt_thal,1);
    ref=[];
    for ct=1:size(BZ.filt_thal{ik,1},1)
        non_nomr=[];epochs_t=[];
        ref=BZ.onset_raw{1,ik}{duration,1};
        env=abs(hilbert(BZ.filt_thal{ik,1}(ct,:)));
        el=500;
        for ii=1:length(ref)
            if ref(ii)>el & ref(ii)+el<length(env)
                epochs_idx(ii,:)=ref(ii)-el:ref(ii)+el;
                epochs_t(ii,:)=(env(ref(ii)-el:ref(ii)+el)-median(env(ref(ii)-el:ref(ii))))./median(env(ref(ii)-el:ref(ii)));
                epochs_ctx(ii,:)=(BZ.env_ctx(ik,(ref(ii)-el:ref(ii)+el))-median(BZ.env_ctx(ik,ref(ii)-el:ref(ii))))./median(BZ.env_ctx(ik,ref(ii)-el:ref(ii)));
            end
        end


        for n=1:(length(env)/100)
            idx_sur=randi([el+1,(length(env)-el)],1,1);
            epochs_idx_sur(n,:)= idx_sur-el:idx_sur+el;
            epochs_t_sur(n,:)= (env(idx_sur-el:idx_sur+el)-median(env(idx_sur-el:idx_sur)))./median(env(idx_sur-el:idx_sur));;
        end
        
        med_thal(ct,:)=median(epochs_t);
        med_sur(ct,:)=median(epochs_t_sur);
        
        clearvars epochs_t epochs_t_sur 
    end
    
    if duration==1
         BZ.pchange_s(ik,:)=mean(med_thal,1);
         BZ.pchange_surr(ik,:)=mean(med_sur,1);
         BZ.pchange_ctxs(ik,:)=mean(epochs_ctx);
    elseif duration==2
        BZ.pchange_l(ik,:)=mean(med_thal,1);
        BZ.pchange_ctxl(ik,:)=mean(epochs_ctx);
    end
        clearvars med_thal med_sur epochs_ctx
end 
end


clearvars -except BZ
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal')
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal')
save 'BZ_opt.mat'

