clearvars -except med_ctx med_thal stats_bursts
thal_long= vertcat(prctile(med_thal,25) ,median(med_thal) ,...
    prctile(med_thal,75));

ctx_long= vertcat(prctile(med_ctx,25) ,median(med_ctx) ,...
    prctile(med_ctx,75));

cd ('C:\Users\wolf5173\Documents\Carolina_code\codes_thal\A1_Thal\mat')
clearvars -except thal_long ctx_long stats_bursts
save SNR_long

thal_sls={[thal_short] ; [thal_long] ; [thal_surr]};
ctx_sl={[ctx_short] ; [ctx_long]};

clearvars -except thal_sls ctx_sl stats_bursts
save SNR