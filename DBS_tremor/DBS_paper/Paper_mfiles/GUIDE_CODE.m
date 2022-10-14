clear; close all
load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','aegean','stone','squash','sapphire','azure');
color_b1=[aegean;stone;blushred];

spiral=0;

if spiral==0
   load('/Users/Carolina/Desktop/oxford/data_code_thesis/DBS/paper_analyses/paper_matfiles/DBS_final_posture_unclust.mat');
else
    load('/Users/Carolina/Desktop/oxford/data_code_thesis/DBS/paper_analyses/paper_matfiles/DBS_final_spiral_unclust.mat');
end


[match_ax]=link_ax(spiral);


%%%%% fx for figures figures

[avg_sup,supp1]=three_by4(out,color_b1);%%% Figure 1A and 2A (select spiral=0 or spiral =1). arcs with significant phases for all patients and axes
  
[p1]=main_psi(out,color_b1); %%%Figure 1B and 2B (select spiral=0 or spiral =1). fig(1) histogram of axes/n trials (Figure 1B and 2B); fig(2) PSI between No Stim main axis and random stim main axis (Figure Supllempentary S 2)

[p1]=axes_count(out); %%% (select spiral=0 or spiral =1) percentage of trials where x y and z are the most predominant

[x,f1,y]=trem_prop(out,match_ax,spiral); %%% coeficcient of variation of unstimualted amp-freq tremor profiles. correlated with deviation of most suppresive effect from level of natural supression 

[avg_sup]=plot_arc(out,match_ax,color_b1); %%% fig(1) arc per patient main axis with sig phases; fig(2) arc per pt all axes

[f1]=modintime_5sec(out,match_ax); %%% modulation in time of most suppresive phase per patient. Change "modu=2" for most amp phase and "axi" for other axes  

%%%other codes

cont_posture.m %%% Figure 3; prolonged phase-specific stimualtion during posture

ACC_metrics.m %%% data for group comparison of ACC spiral data in the 3 stim conditions: NS/HFS/PLS

plotting_spiral_metrics.m %%% means from SPSS (Figure 4.A ; Figure 4.B)

spiral_final.m %%% metrics to make comparisons across quadrants spss & heath plots spiral (Figure 4.C)

power_spiral_posture_3ax.m %%% psd of 3 axes during posture and spiral per patient (not used in paper)

DBS_newmaster.m %%% to generate DBS_final_spiral_unclust.mat and DBS_final_posture_unclust.mat

organise_bins.m %%% spiral data from drawings

freq_filtsig.m  %%% frquency of signal 1 Smr - the signal filtered externally through the digitimer


[dev_freq,sev]=inst_freq_cont_id(s,start,ending,tp_s,tp_e,match_ax,pt); %% p TREMOR SEVERITY PSS
NS_metrics.m %% TREMOR SEVERITY NS
Tseverity_HFS.m %% TREMOR SEVRITY HFS


%%% MORE CODE ON

GUIDE_ADDON.m