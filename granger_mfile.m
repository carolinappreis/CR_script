clear all
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal')
% cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal')
load('BZ_opt.mat');

for duration=1:2
    for ik=1:size(BZ.bua_coh,1);
        ref=[];
        for ct=1:size(BZ.bua_coh{ik,1},1)
            ref=BZ.onset_raw{1,ik}{duration,1};
            clear th; th=BZ.bua_coh{ik,1}(ct,:);
            el=256;
            for ii=1:length(ref)
                if ref(ii)>el & ref(ii):ref(ii)+el<length(th)
                    epochs_idx(ii,:)=ref(ii):ref(ii)+el;
                    epochs_t(ct,ii,:)=th(ref(ii):ref(ii)+el);
                    epochs_ctx(ii,:)=BZ.ctx_coh(ik,(ref(ii):ref(ii)+el));
                end
            end
            
            
            for n=1:(length(ref))
                if ref(ii)>el & ref(ii):ref(ii)+el<length(th)
                idx_sur=randi([el+1,(length(th)-el)],1,1);
                epochs_idx_sur(n,:)= idx_sur:idx_sur+el;
                epochs_t_sur(ct,n,:)= th(idx_sur:idx_sur+el);
                epochs_ctx_sur(n,:)=BZ.ctx_coh(idx_sur:idx_sur+el);
                end
            end
            
            
        end
        
        if duration==1
            BZg.inb_ths{ik,1}=epochs_t;
            BZg.sur_th{ik,1}=epochs_t_sur;
            BZg.sur_ctx{ik,1}=epochs_ctx_sur;
            BZg.inb_ctxs{ik,1}=epochs_ctx;
        elseif duration==2
            BZg.inb_thl{ik,1}=epochs_t;
            BZg.inb_ctxl{ik,1}=epochs_ctx;
        end
%                 if duration==1
%                     BZg.outb_ths{ik,1}=epochs_t;
%                     BZg.outb_ctxs{ik,1}=epochs_ctx;
%                 elseif duration==2
%                     BZg.outb_thl{ik,1}=epochs_t;
%                     BZg.outb_ctxl{ik,1}=epochs_ctx;
%                 end
        clearvars epochs_t epochs_t_sur epochs_ctx epochs_idx epochs_ctx_sur 
    end
end


clearvars -except BZg
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal')
% cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal')
% % save 'BZ_granger.mat'

