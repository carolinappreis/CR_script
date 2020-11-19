
%%% choose between 1) patient 3 (DT) 3) patient 6 (ET) stim at 120 deg 4)
%%% pateint 6 (ET) at 240

% clear
%   load('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_pls_sd/P03_PLS_P')
%   load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_aux_out.mat','rs_mat')
% rs_max=rs_mat(2,1); clear rs_mat; cr=[];

% clear
% load('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_pls_sd/P06_PLS_P1')
%   load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_aux_out.mat','rs_mat')
% rs_max=rs_mat(4,1);cr=[];


clear
 load('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_pls_sd/P06_PLS_P2')
  load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_aux_out.mat','rs_mat') %%% to get the main axis and points of tp_s(tapping start) and tp_e(tapping ending)

rs_max=rs_mat(4,1); cr=1;

%%% -------

[d]=dbs_preprocess(SmrData); samplerateold=d.samplerateold; samplerate=d.samplerate;

co=3; %%% condifition 3 - phase locked

iii=1;
clust=struct; out=struct; start=cell(1,3); ending=cell(1,3); yy=cell(1,3); h_up=cell(1,3); s=struct;

[peak_ax, start, ending, yy, h_up]= dbs_startend(d,samplerateold,samplerate,iii,co,start,ending,yy,h_up);

[s]=dbs_zfiltenv(d,peak_ax,co,iii,s,samplerate);

[b,a]=butter(2,[0.8/(0.5*samplerate) ],'low'); %15
for i=1:3
    s.l_filt{iii,co}(i,:)=filtfilt(b,a,s.raw{iii,co}(i,:));
end


%%% PLS06 P2
if cr==1
ending{1,3}(1,2)=258848;
end


% for i=1:3
%     subplot(3,1,i)
% plot(d.ds_dall(2,:))
% hold on
% plot(s.l_filt{iii,co}(i,:))
% ylim([1.9 2])
% end



y=[sum(s.env{iii,co}(1,1:start{1,3}(1):ending{1,3}(1)));sum(s.env{iii,co}(2,1:start{1,3}(1):ending{1,3}(1)));sum(s.env{iii,co}(3,1:start{1,3}(1):ending{1,3}(1)))];
max_env(1,1)=find(y==(max(y))); clear y
for j=1:size(start{1,3},2)
y=[sum(s.env{iii,co}(1,start{1,3}(j)-3000:start{1,3}(j)));sum(s.env{iii,co}(2,start{1,3}(j)-3000:start{1,3}(j)));sum(s.env{iii,co}(3,start{1,3}(j)-3000:start{1,3}(j)))];
max_env_bf(1,j)=find(y==(max(y))); clear y    
    
    
y=[sum(s.env{iii,co}(1,start{1,3}(j):ending{1,3}(j)));sum(s.env{iii,co}(2,start{1,3}(j):ending{1,3}(j)));sum(s.env{iii,co}(3,start{1,3}(j):ending{1,3}(j)))];
max_env_t(1,j)=find(y==(max(y))); clear y
end


for ax=1:3
figure(1)
plot(1:length(d.ds_dall(2,:)),d.ds_dall(2,:),'Color',[0.5 0.5 0.5])
hold on
plot(s.l_filt{iii,co}(ax,:),'LineWidth',1.5)

figure(2)
plot(1:length(d.ds_dall(2,:)),d.ds_dall(2,:),'Color',[0.5 0.5 0.5])
hold on
plot(1:length(s.filt{iii,co}(ax,:)),s.filt{iii,co}(ax,:),'LineWidth',1.5)
plot(1:length(s.env{iii,co}(ax,:)),s.env{iii,co}(ax,:),'LineWidth',1.5)

for p= 1:length(start{iii,co})
    
%         figure(1)
%     xline(start{iii,co}(p))
%     xline(ending{iii,co}(p))
%     figure(2)
%     xline(start{iii,co}(p))
%     xline(ending{iii,co}(p))
%     
%     
    
    BS(1,p)=mean(s.env_acc{iii,co}(ax,start{iii,co}(p)-3000:start{iii,co}(p)));
    OS(1,p)=mean(s.env_acc{iii,co}(ax,ending{iii,co}(p)-3000:ending{iii,co}(p)));

    for K=1:size(tp_s{p,1},2)
        BHD{p,1}(1,K)=mean(s.env_acc{iii,co}(ax,tp_s{p,1}(K)-3000:tp_s{p,1}(K)));
        AHD{p,1}(1,K)=mean(s.env_acc{iii,co}(ax,tp_e{p,1}(K):tp_e{p,1}(K)+3000));
        figure(1)
        xline(tp_s{p,1}(K),'y.','LineWidth',2)
        xline(tp_e{p,1}(K),'y.','LineWidth',2)
        ylim([1.9 2])
        figure(2)
        xline(tp_s{p,1}(K),'y','LineWidth',0.5)
        xline(tp_e{p,1}(K),'y','LineWidth',0.5)
        xline(tp_s{p,1}(K)-3000,'b','LineWidth',0.5)
        xline(tp_e{p,1}(K)+3000,'b','LineWidth',0.5)
        ylim([-0.05 0.05])
    end
    
%     figure(p+2)
    if K==1
        ax_ch{p,1}(ax,:)=([(BHD{p,1}(1,1)-BS(1,p))./BS(1,p)  (AHD{p,1}(1,1)-BS(1,p))./BS(1,p) (OS(1,p)-BS(1,p))./BS(1,p)])
    else
        ax_ch{p,1}(ax,:)=([(BHD{p,1}(1,1)-BS(1,p))./BS(1,p)  (AHD{p,1}(1,1)-BS(1,p))./BS(1,p) (BHD{p,1}(1,2)-BS(1,p))./BS(1,p) (AHD{p,1}(1,2)-BS(1,p))./BS(1,p) (OS(1,p)-BS(1,p))./BS(1,p)])
        %         xlim([0.5 5.5])
        %         xticks(1:5)
        %         xticklabels({'BS','BHD','HD','AHD','OS'})
    end
end

close all
end



load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','aegean','stone','squash','sapphire','azure','ax2','ax3');

color_b1=[stone;ax2;ax3];

for i=1:length(start{iii,co})
f1=figure(i);
set(f1,'color','w')
b=bar(ax_ch{i,1}','FaceAlpha',[0.75],'EdgeColor','none')
b(1).FaceColor=stone;
b(2).FaceColor=ax2;
b(3).FaceColor=ax3;
box('off')
ylim([-1 1])
yticks(-1:0.2:1)
ylabel({'Change in tremor severity'})
legend('tracked axis','axis2','axis3')
legend('Orientation','horizontal')
legend('Location','southwest')
legend('boxoff')
set(gca,'xtick',[]);
set(gca,'xcolor',[1 1 1])
set(gca,'FontSize',14);
end


for i=1:3
f1=figure(i);
set(f1,'color','w')
plot(1:5,ax_ch{i,1}(1,:),'--','Color',stone,'LineWidth',1)
hold on
stem(1:5,ax_ch{i,1}(1,:),'.', 'LineWidth',4,'MarkerSize',20,'Color',stone)
xlim([-0.5 6.5])

%%% ib case of just one tap use

% % plot(1:3,ax_ch{i,1}(1,:),'--','Color',stone,'LineWidth',1)
% % hold on
% % stem(1:3,ax_ch{i,1}(1,:),'.', 'LineWidth',4,'MarkerSize',20,'Color',stone)
% % xlim([-0.5 4.5])

box('off')
ylim([-1 1])
yticks(-1:0.2:1)
ylabel({'Change in tremor severity'})
set(gca,'xtick',[]);
set(gca,'xcolor',[1 1 1])
set(gca,'FontSize',14);
end






%%% find mhd mhu points
% take mean evelope second before hand up before mhd during handdown second
% after mhup

% plot(1:length(d.ds_dall(2,:)),d.ds_dall(2,:))
% hold on
% plot(1:length(d.ds_dall(1,:)),d.ds_dall(1,:))

%%%% points tapping saved in mat files

% % % % %%%03 pls

% % % % s=[];
% % % % e=[];

% % % % %%%06 pls p1
%%% BEGINING AND END STIM AT TAPPING
% % tp_s=({[50490 91804];[236459 264248];[434366]});
% % tp_e=({[51970 94613]; [237689 266437];[436481]});

%%% MATCHED TO DC SIGNAL
% % tp_s=({[44168 85467];[233126 261009];[428515]});
% % tp_e=({[55987 98462]; [240524 267991];[438380]});


% % % % %%%06 pls p2
%%% BEGINING AND END STIM AT TAPPING
% %  tp_s=({[71918 91887];[220416 261695]});
% %  tp_e=({[73931 93599];[223301 263252]});

%%% MATCHED TO DC SIGNAL
% % %  tp_s=({[65517 89529];[217603]});
% % % tp_e=({[76628 96325];[224831]});

%%%% p03 PLS 
%%% MATCHED TO DC SIGNAL
% % % tp_s=({[111279 143311];[339162 368494];[557432 585324]});
% % % tp_e=({[116458 148097];[343612 373744];[561820 589888]});



% for j=1:size(start{1,3},2)
% x=[(s.filt{iii,co}(1,start{1,3}(j):ending{1,3}(j)));(s.filt{iii,co}(2,start{1,3}(j):ending{1,3}(j)));(s.filt{iii,co}(3,start{1,3}(j):ending{1,3}(j)))];
% [pc, score, latent, tsquare, explained] = pca(x');
%  max_pc(j,:)=find(pc(:,1)==max(pc(:,1))); clear x
% end
