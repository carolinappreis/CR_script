

%%% choose between 1) patient 3 (DT) 3) patient 6 (ET) stim at 120 deg 4)
%%% pateint 6 (ET) at 240

% clear
%   load('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_pls_sd/P03_PLS_P')
%   load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_aux_out.mat','rs_mat')
% rs_max=rs_mat(2,1); clear rs_mat; cr=[];

clear
load('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_pls_sd/P06_PLS_P1')
  load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_aux_out.mat','rs_mat')
rs_max=rs_mat(4,1);cr=[];


% clear
%  load('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_pls_sd/P06_PLS_P2')
%   load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_aux_out.mat','rs_mat') %%% to get the main axis and points of tp_s(tapping start) and tp_e(tapping ending)

rs_max=rs_mat(4,1); cr=1;

%%% -------

[d]=dbs_preprocess(SmrData); samplerateold=d.samplerateold; samplerate=d.samplerate;

co=3; %%% condifition 3 - phase locked

iii=1;
clust=struct; out=struct; start=cell(1,3); ending=cell(1,3); yy=cell(1,3); h_up=cell(1,3); s=struct;

[peak_ax, start, ending, yy, h_up]= dbs_startend(d,samplerateold,samplerate,iii,co,start,ending,yy,h_up);

[s]=dbs_zfiltenv(d,peak_ax,co,iii,s,samplerate);

dum=[];
for g=1:size(tp_s,1)
    for m=1:size(tp_s{g,1},2)
       dum=[dum tp_s{g,1}(1,m):tp_e{g,1}(1,m)];
    end
end
hands=h_up{1,3}(1,:);
[gh,rh]=intersect(hands,dum);
hands(1,rh)=NaN; hands1=hands(~isnan(hands)); hands1=sort(hands1,'ascend');

% t=1:length(s.env{1,3}(1,:));
% plot(t,s.env{1,3}(1, :))
% hold on
% plot(t(1,h_up{1,3}),s.env{1,3}(1,h_up{1,3}))
% plot(t(1,dum),s.env{1,3}(1,dum))
% plot(t(hands1),s.env{1,3}(1,hands1),'--')

all_idx=[1:start{1,3}(1) hands1 ending{1,3}(end):length(s.env{1,3})];

% plot(t,s.env{1,3}(1,:))
% hold on
% plot(t(1,all_idx),s.env{1,3}(1,all_idx))

envelope=s.env{1,3}(1,:);

env_s=s.env{1,3}(1,1:start{1,3}(1))

unstable=find(envelope(1:start{1,3}(1))<(mean(envelope(all_idx))-std(envelope(all_idx))));