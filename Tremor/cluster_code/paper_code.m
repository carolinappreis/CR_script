%%%%%%%%%%%%%%%%%%% MEDIAN NERVE PHASE-SPECIFIC STIM CODE %%%%%%%%%%%%%%%%%%

%RESULTS PAPER
cd('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data')
cd('/Users/Carolina/Documents/GitHub/CR_script/Tremor/cluster_code')
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/cluster_out_mc.mat');

        %%% 1) plots & stats stim ARCs and mediansplit stim ARCs with clustering 
        [r,p]=plots_c(out);


        %%% 2) plots & stats tremor power during non-stim / superssive stim effects and amplifying stim effects 
        [out]=psd_rs_indseg(out); 
        % avg pwelch of 5 sec segments amp and sup - 
                %%% pwelch from 5 sec seg of 50 000 non-stim surrogate comes from
                %%% ppx_nseg.m

         %%% 3) plots & stats group stim effects vs non stim effects      
        [p_amp,p_sup]=group_arcs(out);

         %%% 4) plots & stats group stim effects vs non stim effects after median
         %%% split
        run('median_split_group.m');


%%%% 1)-3) need variable 'out' from:
run('cluster_master_clean.m') 
    %%includes functions:
    % [d]=preprocess(SmrData);
    % [h]=twitch(iii,d,samplerateold,samplerate,co);
    % [afilt, bfilt, start, ending, yy, h_up]= startend(d,samplerateold,samplerate,iii,co,start,ending,yy,h_up);
    % [s]=zfiltenv(d,bfilt,afilt,co,iii,s);
    % **[clust,out]=clustering2(out,iii,clust,start,ending,yy);
    % ** [out]=mod_c(clust,out,co,iii,s,cohort);
    % **=will ask for output from code run('aux_cluster_master.m') named aux_out.mat
