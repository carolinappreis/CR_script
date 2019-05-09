clear all
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal')
load('BZ.mat')



for ik=1:size(BZ.env_ctx,1);
 env=(BZ.env_ctx(ik,:)-median(BZ.env_ctx(ik,:)))./median(BZ.env_ctx(ik,:));
 BZ.Nonset_raw{1,ik}=bursts(env);
 BZ.Noffset_raw{1,ik}=bursts_off(env);
 Ecogfiltered=BZ.filt_ctx(ik,:);
 BZ.Nonset_phase_al{1,ik}=bursts_aligned(env,Ecogfiltered);
 BZ.Noffset_phase_al{1,ik}=bursts_aligned_off(env,Ecogfiltered);
end




    

    

