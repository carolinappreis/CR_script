close all
clear all


load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','aegean','stone','squash','sapphire','azure','ax2','ax3');

color_b1=[stone;ax2;ax3];

spiral=0;[match_ax]=link_ax(spiral);

%%% choose between 1) patient 3 (DT) 3) patient 6 (ET) stim at 120 deg 4)
%%% pateint 6 (ET) at 240
% 
% load('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_pls_sd/P03_PLS_P')
% cr=[];pt=2;

load('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_pls_sd/P06_PLS_P1')
cr=[];pt=4
% 
% load('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_pls_sd/P06_PLS_P2')
% cr=1;pt=4

%%% -------

[d]=dbs_preprocess(SmrData); samplerateold=d.samplerateold; samplerate=d.samplerate;

co=3; %%% condifition 3 - phase locked

iii=1;
out=struct; start=cell(1,3); ending=cell(1,3); yy=cell(1,3); h_up=cell(1,3); s=struct;

[peak_ax, start, ending, yy, h_up]= srtend(d,samplerateold,samplerate,iii,co,start,ending,yy,h_up,spiral);

[s]=zfiltenv_simple(d,peak_ax,co,iii,s,samplerate);

[f1]=in_time(iii,s,co,samplerate,start,ending,tp_s,tp_e,color_b1,cr,match_ax,pt);

% [f1]=pc_ax(s,tp_s,tp_e,cr,start,ending,color_b1,match_ax,pt);