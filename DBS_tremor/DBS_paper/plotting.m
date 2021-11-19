clear; close all
load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','aegean','stone','squash','sapphire','azure');
color_b1=[aegean;stone;blushred];

spiral=1;

if spiral==0
   load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_final_posture.mat');
else
    load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_final_spiral.mat');
end
[match_ax]=link_ax(spiral);




%%%%% figures
% [avg_sup]=plot_arc(out,match_ax,color_b1);
%
% [f1]=sup_envtime(out,match_ax,spiral);
% [avg_sup]=three_by4(out,color_b1);
% %   
%  [p1]=main_psi(out,color_b1);
%
[f1,y]=trem_prop(out,match_ax,spiral);
