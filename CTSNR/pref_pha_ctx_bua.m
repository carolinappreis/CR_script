function [fig]=pref_pha_ctx_bua(coh_filts,name)

% r=0;

for pr=1:size(coh_filts,1)
    ctx=coh_filts{pr,1}(1,:);
    thal=coh_filts{pr,1}(2:end,:);
    
    env=abs(hilbert(ctx));
    %         env=smoothdata(env,'movmean',50);
    %         [onset1,offset1]=bursts(env);
    
    [onset,offset]=bursts_aligned(env,ctx);
    onset1=onset; clear onset
    offset1=offset; clear offset
    onset=horzcat(onset1{:});
    offset=horzcat(offset1{:});
    
    ctx_phase=angle(hilbert(ctx));
    all_ct=[];
    for ct=1:size(thal,1)
        clear thal_phase non_norm
        thal_phase=angle(hilbert(thal(ct,:)));
        non_norm=wrapToPi(ctx_phase-thal_phase);
        
        clear epochs_t
        for jj=1:length(onset)
            el=500;
            if onset(jj)>el && onset(jj)+el<length(non_norm)
                epochs_t(1,jj)=circ_mean((non_norm(onset(jj)-el:onset(jj)+el))');
            end
            
        end
        
        all_ct=[all_ct epochs_t];
        
        %             clear epochs_idx_sur epochs_t_sur
        %             for n=1:(length(non_norm)/1000)
        %                 idx_sur=randi([el+1,(length(non_norm)-el)],1,1);
        %                 epochs_t_surr(pr,ct,n)=circ_mean((non_norm(idx_sur-el:idx_sur+el))'); clear idx_sur
        %             end
        
    end
    phase_m(pr,:)=circ_mean(all_ct'); 
    phase_std(pr,:)=circ_std(all_ct'); clear all_ct
    phall(pr,1)=circ_mean(non_norm');
    phall(pr,2)=circ_r(non_norm');
    phall(pr,3)=(circ_rtest(non_norm')==0);
     
    
    clearvars -except pr coh_filts el name r phall  phase_m phase_std
end
in_b=360+rad2deg(circ_mean(phase_m))
all=360+rad2deg(circ_mean(phall(:,1)))
end