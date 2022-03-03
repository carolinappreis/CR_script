close all
clear all


load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','aegean','stone','squash','sapphire','azure','ax2','ax3');
color_b1=[stone;ax2;ax3];

spiral=0;[match_ax]=link_ax(spiral);

%%% choose between:

%%% 1) patient 3 (DT) 
load('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_pls_sd/P03_PLS_P')
cr=[];pt=2;

%%%  2) patient 6 (ET) stim at 120 deg (most sup)
% % % load('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_pls_sd/P06_PLS_P1')
% % % cr=[];pt=4;

%%%  3)pateint 6 (ET) stim at 240 deg (2nd most sup)
% load('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_pls_sd/P06_PLS_P2')
% cr=1;pt=4


%%% used in paper: 1) and 2)

%%% -------

[d]=dbs_preprocess(SmrData); samplerateold=d.samplerateold; samplerate=d.samplerate;

co=3; %%% condifition 3 - phase locked

iii=1;
out=struct; start=cell(1,3); ending=cell(1,3); yy=cell(1,3); h_up=cell(1,3); s=struct;

[peak_ax, start, ending, yy, h_up]= srtend(d,samplerateold,samplerate,iii,co,start,ending,yy,h_up,spiral);

[s]=zfiltenv_simple(d,peak_ax,co,iii,s,samplerate);


[f1]=in_time(iii,s,co,samplerate,start,ending,tp_s,tp_e,color_b1,cr,match_ax,pt,d);

[total]=coef_share(s,start,ending,tp_s,tp_e,match_ax,pt);
% [f1]=pc_ax(s,tp_s,tp_e,cr,start,ending,color_b1,match_ax,pt);


% % phase=s.phase{iii,3};
% % phase_difs=[wrapToPi(phase(1,:)-phase(2,:));wrapToPi(phase(1,:)-phase(3,:));wrapToPi(phase(2,:)-phase(3,:))];
% % b=start{iii,3};
% % e=ending{iii,3};
% % for ax=1:2
% %     for st=1:length(b)
% %         m(ax,st)=circ_r((phase_difs(ax,b(st):e(st)))');
% %     end
% % end
% % 
% % [mean(m') ; std(m')]


