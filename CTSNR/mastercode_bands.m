clear all
close all
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/final_mats')
load('BZ_bua.mat');

%filtering probe signals and cortex in the beta band (+-5Hz peak coherence) if the summed coherence
%between 15-35Hz is more than 10% of the coherence in all freqs.

samprate=1000;
for pr=1:size(bua,1)
    m=2;
    for ct=2:size(bua{pr,1},1)
        [Pxx_ind,F_ind]=mscohere(bua{pr,1}(1,:),bua{pr,1}(ct,:),samprate,[],samprate,samprate);
        Pxx_ind_beta=Pxx_ind(16:36);filtrange=14+find(Pxx_ind_beta==max(Pxx_ind_beta));
        [b,a]=butter(2,[(filtrange-5)/(0.5*samprate) (filtrange+5)/(0.5*samprate)],'bandpass');
        if (sum(Pxx_ind(16:36))/sum(Pxx_ind(1:end)))>0.1
            % band=[35 60];
             band=[0.5 3];
%              band=[4 12];
            [c,d]=butter(2,[(band(1))/(0.5*samprate) (band(end))/(0.5*samprate)],'bandpass');
            coh_filts{pr,1}(1,1:length(bua{pr,1}(1,:))) = filtfilt(b,a,bua{pr,1}(1,:));        
            coh_filts{pr,1}(m,1:length(bua{pr,1}(1,:))) = filtfilt(c,d,bua{pr,1}(ct,:));
            m=m+1;
        end
        clear Pxx_ind_beta filtrange
    end
    ctx(pr,:)=filtfilt(b,a,bua{pr,1}(1,:));
end
coh_filts=coh_filts(~cellfun('isempty',coh_filts));
ctx=ctx((~cellfun('isempty',coh_filts)),:);

% [env_var]=burst_var(ctx);
%%%%%%%%

clearvars -except coh_filts name


% CHANGE IN BETA AMPLITUDE
% [fig]=ctx_to_sub(coh_filts,name);
% [fig]=sub_to_ctx(coh_filts,name);

% CHANGE IN BETA PHASE COUPLING
% [fig]=sub_to_ctx_psi(coh_filts,name);
% [fig]=ctx_to_sub_psi(coh_filts,name);


% PHASE SLIPS IN ALIGNMENT AND INSTANTANEOUS FREQUENCY
% [fig]=pha_slip_ctx_subctx(coh_filts,name)
% [fig]=pha_slip_ctx(coh_filts,name)
% [fig]=pha_slip_sub(coh_filts,name)

% [fig]=pha_slip_subctx_ctx(coh_filts,name)
% [fig]=pha_slip_ctx_subref(coh_filts,name)
% [fig]=pha_slip_sub_subref(coh_filts,name)

% BURSTS OVERLAP
% [fig]=subctx_ovl(coh_filts,name);
% [fig]=ecog_sub_ovl(coh_filts,name);
% [fig]=sub_ecog_ovl(coh_filts,name);

