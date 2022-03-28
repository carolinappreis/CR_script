%%%%%%%%%%%%%%%%%%% MEDIAN NERVE PHASE-SPECIFIC STIM CODE %%%%%%%%%%%%%%%%%%
%%% .m files can be found in my Github repositrory: https://github.com/carolinappreis/CR_script.git
%%% under the folder tremor > cluster_code

%RESULTS PAPER
        %%% 1) plots & stats stim ARCs and mediansplit stim ARCs with
        %%% clustering (fig.1paper = fig 1code; fig.3A= fig 2code; fig3B= fig 3code)
        [r,p]=plots_c(out); 
        %%% 2) plots & stats group stim effects vs non stim effects      
        [p_amp,p_sup]=group_arcs(out);

         %%% 3) plots & stats tremor power during non-stim / superssive stim effects and amplifying stim effects 
        [out]=psd_rs_indseg(out); % asks for pxx_ns_indivseg.mat , ehich comes from  'ppx_nseg.m'
            % avg pwelch of 5 sec segments amp and sup - pwelch from 5 sec seg
            % of 50 000 non-stim surrogate 

   

%%%% 1) to 3) need variable 'out' found in 'cluster_out_mc.mat' which was saved from:
run('cluster_master_clean.m') 
    %%includes functions:
    % [d]=preprocess(SmrData);
    % [h]=twitch(iii,d,samplerateold,samplerate,co);
    % [afilt, bfilt, start, ending, yy, h_up]= startend(d,samplerateold,samplerate,iii,co,start,ending,yy,h_up);
    % [s]=zfiltenv(d,bfilt,afilt,co,iii,s);
    % **[clust,out]=clustering2(out,iii,clust,start,ending,yy);
    % ** [out]=mod_c(clust,out,co,iii,s,cohort);
    % **=will ask for output from code run('aux_cluster_master.m') named aux_out.mat
