
clear all; close all
load('/Users/Carolina/Documents/GitHub/CR_script/DBS_tremor/DBS_paper/all_cond_spiral.mat')
load('/Users/Carolina/Documents/GitHub/CR_script/DBS_tremor/DBS_paper/comb_pct.mat')


for i=1:3
    bar_tremor(i,:)=avg_t_sev(i,:,1)./(pct(i,2).*(10*9.81/0.5));
    
    signal=ufsignal_all{i,1}; samplerate =1000;
    
    [Pxx,F] = pwelch(signal, samplerate, [], round(samplerate), samplerate);
    range = Pxx(3:10);
    curve(iii-1,cond,:)=Pxx./(sum(Pxx(100:500)));
    peak(i,:)=max(squeeze(curve(iii-1,cond,3:10)));
    
end

figure
bar(bar_tremor)

figure
for ii=1:3
    for cc=4
        subplot(1,3,ii)
        plot(F(3:11),(squeeze(Pcurve(ii,cc,3:11)))./ peak(ii))
        hold on
    end
    point(ii,:)=Ppeak(ii,:)./peak(ii);
    for i=2:4
        change(ii,i-1,:)=round(((point(ii,i,:)-point(ii,1,:))./point(ii,1,:)*100),0);
    end
    box('off')
    xlim([2 10])
end

figure
bar(change)
box('off')
ylabel('% supression from baseline')
xlabel('patients')



load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','aegean','stone','squash','sapphire','azure');
color_b1=[aegean;stone;blushred];

avg_dur=(squeeze(nanmean(dur_b,1))); avg_int=(squeeze(nanmean(int_b,1))); avg_nr=mean(nr_above,1)
avg_dur(1:3,:)=NaN;
figure()
bar(avg_dur([1 3 4],1))
hold on
errorbar([1:3],avg_dur([1 3 4],1),avg_dur([1 3 4],2),'.','Color',color_b1(1,:),'LineWidth',2)
box('off')
figure()
bar(avg_int([1 3 4],1))
hold on
errorbar([1:3],avg_int([1 3 4],1),avg_int([1 3 4],2),'.','Color',color_b1(1,:),'LineWidth',2)
box('off')


data=dur_b;
f1=figure()
for iii=1:3

    subplot(1,3,iii)
    err=squeeze(data(iii,:,2));
    x=squeeze(data(iii,:,1));
    bar([x],'EdgeColor','none','FaceColor',color_b1(1,:),'FaceAlpha',0.5)
    hold on
    %     errorbar([1:4],x,err,'.','Color',color_b1(1,:),'LineWidth',2)
    % ylim([0 500])
    % yticks(0:100:500)
    xticklabels({'NS','RS','PLS','HFS'})
    if data==dur_b
        ylabel('burst duration')
    else
        ylabel('int burst interval')
    end
    box('off')
    set(gca,'FontSize',14)
    clear x err
end
f1.OuterPosition= [1,100,1000,300];
set(f1,'color','w');



% for i=1:3
% psi_pls(i,:)=[mean(squeeze(psi(i,3,:))) std(squeeze(psi(i,3,:)))]
% end
% for i=1:3
% subplot(1,3,i)
% bar(squeeze(psi(i,:,:)))
% end
