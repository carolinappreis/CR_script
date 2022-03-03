clear all
close all
cd('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA')
cd('/Users/Carolina/Documents/GitHub/CR_script/DBS_tremor/DBS_paper')
clear; close
cohort = [ 1 3 4 6];
cond={'NS';'RS'};
clust=struct; out=struct; start=cell(10,1); ending=cell(10,1); yy=cell(10,3); h_up=cell(10,3); s=struct; freq_bl=[]; amp_bl=[];
rng('default')
gen=(rng);
spiral=0;

if spiral==0
    load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_final_posture.mat');
else
    load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_final_spiral.mat');
end

co=1;

for iii=1:4
load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA/P0',num2str(cohort(iii)),'_',num2str(cond{co,1}),'.mat'))
[d]=dbs_preprocess(SmrData); samplerateold=d.samplerateold; samplerate=d.samplerate;
[peak_ax, start, ending, yy, h_up]= srtend(d,samplerateold,samplerate,iii,co,start,ending,yy,h_up,spiral);
[s]=zfiltenv_simple(d,peak_ax,co,iii,s,samplerate);
phase=s.phase{iii,1};
phase_difs=[wrapToPi(phase(1,:)-phase(2,:));wrapToPi(phase(1,:)-phase(3,:));wrapToPi(phase(2,:)-phase(3,:))];
b=start{iii,1};
e=ending{iii,1};
for ax=1:2
    for st=1:length(b)
        m(ax,st)=circ_r((phase_difs(ax,b(st):e(st)))');
    end
end
load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','aegean','stone','squash','sapphire','azure');
color_b1=[aegean;stone;blushred];
p1=figure(1)
subplot(1,4,iii)
bar(mean(m'),'EdgeColor','none','FaceColor',color_b1(1,:),'FaceAlpha',0.5)
hold on
ylim([0 1])
yticks(0:0.25:1)
xticklabels({['Z,Y'],['Z,X']})
xlabel ('Axes')
ylabel('PSI')
box('off')
set(gca,'FontSize',14)
end
p1.OuterPosition= [1,100,1000,300];
set(p1,'color','w');