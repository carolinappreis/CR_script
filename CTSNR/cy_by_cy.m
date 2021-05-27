
function [fig]=cy_by_cy(data,ecog,name)

srn=1000;
tt=0;
all=[];
allall=[];
for j=1:size(data,1)
    
    clearvars -except name ecog data srn j b a  check vec_lg pref_pha tt all allall
    ctx=ecog(j,:);
    Ecogfiltered=ecog(j,:);
    env=abs(hilbert(Ecogfiltered));
    data_ones=find(data(j,:)==1);
    hp=wrapToPi(angle(hilbert(Ecogfiltered)));
    ang=hp(data_ones);
    
    
    cy_bursts=cycles_10(env,Ecogfiltered);
    block=cell2mat(cy_bursts);
    for d1=1:size(block,1)
        for d2=1:size(block,2)
            if d2+1<length(block(d1,:))
                epoch=block(d1,d2):block(d1,d2+1);
                l=find(data(j,epoch)==1);
                %avg spikes
                if ~isempty(epoch(l))
                    if (length(epoch))>1
                        pha_b(d1,d2)=circ_mean(hp(epoch(l))'); %%% picking just the first spike in a cycle l(1)
                        all=[all circ_mean(hp(epoch(l))')];
                    else
                        pha_b(d1,d2)=hp(epoch(l(1))); %%% picking just the first spike in a cycle l(1)
                        all=[all hp(epoch(l(1)))];
                    end
                    allall=[allall hp(epoch(l))];
                else
                    pha_b(d1,d2)=NaN;
                end
                
                % % %                     % last spike for those with more than one
                % % %                     if ~isempty(epoch(l))
                % % %                         if (length(epoch))>1
                % % %                         pha_b(d1,d2)=hp(epoch(l(end))); %%% picking just the first spike in a cycle l(1)
                % % %                         all=[all hp(epoch(l(end)))];
                % % %                         else
                % % %                         pha_b(d1,d2)=hp(epoch(l(1))); %%% picking just the first spike in a cycle l(1)
                % % %                         all=[all hp(epoch(l(1)))];
                % % %                         end
                % % %                         allall=[allall hp(epoch(l))];
                % % %                     else
                % % %                         pha_b(d1,d2)=NaN;
                % % %                     end
                % % %
                % % %                     % first spike all
                % % %                     if ~isempty(epoch(l))
                % % %                         pha_b(d1,d2)=hp(epoch(l(1))); %%% picking just the first spike in a cycle l(1)
                % % %                         all=[all hp(epoch(l(1)))];
                % % %                         allall=[allall hp(epoch(l))];
                % % %                     else
                % % %                         pha_b(d1,d2)=NaN;
                % % %                     end
            end
        end
    end
    
    % %     n=40; (bu(1:10,1) if we want fixed number of burst with spikes per
    % cycle
    for x=1:size(pha_b,2)
        bu=pha_b(:,x); bu=bu(~isnan(bu));
        check1(1,x)=length(bu);
        vl(1,x)=circ_r(bu);
        pp(1,x)=circ_mean(bu); clear bu bu1
    end
    
    %     if min(check1)>size(block,1)/4
    if min(check1)>=5
        tt=tt+1;
        vec_lg(tt,:)=vl;
        pref_pha(tt,:)=pp;
        check(tt,:)=check1;
        
% %             figure(1);
% %             polarplot([0 circ_mean(ang')], [0, circ_r(ang')],'linewidth',2)
% %             hold on
% % %           polarhistogram(ang,'BinWidth',2*pi/12)
% %             idx(tt)=j;
% %             figure(2)
% %             polarplot([0 circ_mean(pp')], [0, mean(vl)],'linewidth',2)
% %             hold on
    end
    
    
end

figure(2)
hold on
polarplot([0 circ_mean((circ_mean(pref_pha))')], [0,mean(mean(vec_lg)) ],'r','linewidth',6)

[~,~,stats]=kruskalwallis(vec_lg);
[c,m,h,nms]=multcompare(stats,'CType','dunn-sidak');


vec_m=nanmean(vec_lg,1);
ang_m=circ_mean(pref_pha); rad2deg(circ_mean(ang_m'))
ang_std=circ_std(pref_pha); rad2deg(circ_mean(ang_std'))
ang_sem= (circ_mean(pref_pha)+circ_std(pref_pha))./sqrt(size(vec_lg,1));rad2deg(circ_mean(ang_sem'))
% polarhistogram(pref_pha,'BinWidth',2*pi/12)

%% new figures 5 cycles

load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','squash','blood','sky','aegean','seafoam');
if name=='bz'
    color_b={squash blood};
else
    color_b={sky aegean};
end



fig=figure()
cy_plots=1:11; % first 5
for i=1:length(cy_plots)
    if ~isempty(intersect(cy_plots(i), [1:5]))
        cl=color_b{1,1};
    elseif cy_plots(i)==6
        cl=seafoam;
    else
        cl=color_b{1,2};
    end
    if i==1
        polarplot([0 ang_m(1,cy_plots(i))], [0, vec_m(1,cy_plots(i))],'Color',cl,'linewidth',1,'MarkerIndices',[3],'Marker','d')
    else
        polarplot([ang_m(1,cy_plots(i-1)) ang_m(1,cy_plots(i))], [vec_m(1,cy_plots(i-1)), vec_m(1,cy_plots(i))],'Color',cl,'linewidth',1,'MarkerIndices',[2],'Marker','d')
    end
    hold on
end
legend({'-5cy';'-4cy';'-3cy' ;'-2cy';'-1cy';'burst onset';'+1cy';'+2cy';'+3cy';'+4cy';'+5cy'})
legend('boxoff')
box('off')
set(gca,'FontSize',12)
fig.Units = 'centimeters';
fig.OuterPosition= [10, 10, 15, 15];
fig.Color='w';





fig=figure()
% plot(nanmean(vec_lg,1),'-d','LineWidth',1.5,'Color',color_b{1,2})
x=1:11;
y=mean(vec_lg);
err=std(vec_lg);
% err=(mean(vec_lg)+std(vec_lg)./sqrt(size(vec_lg,1)));
errorbar(x,y,err,'-d','LineWidth',1.5,'Color',color_b{1,2})
hold on
xline(6,'--',{'burst onset'},'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','Color','r','LineWidth',2)
xlim([0 11])
ylim([0 0.3])
yticks(0:0.1:0.3)
xlabel('Number of {\beta} cycles')
xticks([1:1:11])
xticklabels ({'-5','-4','-3','-2','-1','0','1','2','3','4','5'})
ylabel('Vector length')
box('off')
set(gca,'FontSize',12)
fig.Units = 'centimeters';
fig.OuterPosition= [10, 10, 10, 10];
fig.Color='w';

%% old figures 10 cycles
fig=figure()
cy_plots=6:16; % first 5
for i=1:length(cy_plots)
    if ~isempty(intersect(cy_plots(i), [1:10]))
        cl=color_b{1,1};
    elseif cy_plots(i)==11
        cl=seafoam;
    else
        cl=color_b{1,2};
    end
    if i==1
        polarplot([0 ang_m(1,cy_plots(i))], [0, vec_m(1,cy_plots(i))],'Color',cl,'linewidth',1,'MarkerIndices',[3],'Marker','d')
    else
        polarplot([ang_m(1,cy_plots(i-1)) ang_m(1,cy_plots(i))], [vec_m(1,cy_plots(i-1)), vec_m(1,cy_plots(i))],'Color',cl,'linewidth',1,'MarkerIndices',[2],'Marker','d')
    end
    hold on
end
legend({'-5cy';'-4cy';'-3cy' ;'-2cy';'-1cy';'burst onset';'+1cy';'+2cy';'+3cy';'+4cy';'+5cy'})
legend('boxoff')
box('off')
set(gca,'FontSize',12)
fig.Units = 'centimeters';
fig.OuterPosition= [10, 10, 15, 15];
fig.Color='w';



% % fig=figure()
% % % plot(nanmean(vec_lg,1),'-d','LineWidth',1.5,'Color',color_b{1,2})
% % plot(mean(vec_lg),'-d','LineWidth',1.5,'Color',color_b{1,2})
% % hold on
% % xline(11,'--',{'burst onset'},'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','Color',[0.5 0.5 0.5],'LineWidth',2)
% % xlim([0 22])
% % xlabel('Number of {\beta} cycles')
% % xticks([1:2:21])
% % xticklabels ({'-10','-8','-6','-4','-2','0','2','4','6','8','10'})
% % ylabel('Vector length')
% % box('off')
% % set(gca,'FontSize',12)
% % fig.Units = 'centimeters';
% % fig.OuterPosition= [10, 10, 10, 10];
% % fig.Color='w';


fig=figure()
% plot(nanmean(vec_lg,1),'-d','LineWidth',1.5,'Color',color_b{1,2})
x=1:21;
y=mean(vec_lg);
err=std(vec_lg);
% err=(mean(vec_lg)+std(vec_lg)./sqrt(size(vec_lg,1)));
errorbar(x,y,err,'-d','LineWidth',1.5,'Color',color_b{1,2})
hold on
xline(11,'--',{'burst onset'},'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','Color',[0.5 0.5 0.5],'LineWidth',2)
xlim([0 22])
ylim([0 0.6])
yticks(0:0.2:0.6)
xlabel('Number of {\beta} cycles')
xticks([1:2:21])
xticklabels ({'-10','-8','-6','-4','-2','0','2','4','6','8','10'})
ylabel('Vector length')
box('off')
set(gca,'FontSize',12)
fig.Units = 'centimeters';
fig.OuterPosition= [10, 10, 10, 10];
fig.Color='w';

%%% cycle eleven is burst onset
% fig=figure()
% cy_plots=7:15;
% for i=1:length(cy_plots)
%     subplot(1,length(cy_plots),i,polaraxes)
%     for un=1:size(vec_lg,1)
%         polarplot([0 pref_pha(un,cy_plots(i))], [0, vec_lg(un,cy_plots(i))],'linewidth',2)
%         rlim([0 0.5])
%         hold on
%     end
% end
%
% fig.Units = 'centimeters';
% fig.OuterPosition= [10, 10, 40, 10];
% fig.Color='w';


end

% % %  for single units. fig.1 - compare pref angle all recording versus 1st spike
% % % per cycle; fig.2 - compare pref angle all recording versus all spikes
% % % per cycle
% % figure(1)
% % subplot(1,2,1)
% % polarhistogram(ang,'BinWidth',2*pi/12)
% % hold on
% % polarhistogram(all,'BinWidth',2*pi/12)
% % subplot(1,2,2)
% % polarplot([0 circ_mean(ang')], [0,circ_r(ang') ],'linewidth',6)
% % hold on
% % polarplot([0 circ_mean(all')], [0,circ_r(all') ],'linewidth',6)
% % 
% % figure(2)
% % subplot(1,2,1)
% % polarhistogram(ang,'BinWidth',2*pi/12)
% % hold on
% % polarhistogram(allall,'BinWidth',2*pi/12)
% % subplot(1,2,2)
% % polarplot([0 circ_mean(ang')], [0,circ_r(ang') ],'linewidth',6)
% % hold on
% % polarplot([0 circ_mean(allall')], [0,circ_r(allall') ],'linewidth',6)