%%%% GUIDE FOR ADD_ON CODES



%%%%%    POSTURE

%%% ------ main tremor axis:

%%% random stim: 

[p1]=main_psi(out,color_b1).m %%% run GUI_CODE.m and select spiral==0; fig.1


%%% continous:

cont_post.m %%% run until [s]=zfiltenv_simple(d,peak_ax,co,iii,s,samplerate);
[total]=twomet_coef_share(s,start,ending,tp_s,tp_e,match_ax,pt); %%% choose if mesure is pca or sum env


%%% ------ inst freq:

%%% random stim:

inst_freq_plots_rs.m 


%%% continous:

[dev_freq,sev]=inst_freq_cont_id(s,start,ending,tp_s,tp_e,match_ax,pt); %% per patient ---- FREQ AND TREMOR SEVERITY
[dev_freq]=inst_freq_cont(s,start,ending,tp_s,tp_e,match_ax,pt); 



%%%%%    SPIRAL
%%% ------ main tremor axis:

%%% random stim: 

[p1]=main_psi(out,color_b1).m %%% %%% run GUI_CODE.m and select spiral==1; fig.1


%%% continous: ----- spatial variability

ACC_metrics.m %%% and  %% plotting_spiral_metrics
organise_bins.m %%% to plot spiral with tremor power 

%%% ------ inst freq:

%%% random stim:

inst_freq_plots_rs.m 


%%% continous:

ACC_instfreq_spiral.m
bin_data_freq.m %%% to plot spirals with tremor freq deviation and uncomment, across quadrants. 