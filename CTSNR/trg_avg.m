function [fig]=trg_avg(data,ecog,name)

srn=1000;
[b,a]=butter(2,[15/(0.5*srn) 35/(0.5*srn)],'bandpass');
tt=0;
for j=1:size(data,1)
    
    clearvars -except j b a ecog data srn rec_pa rec_npa count_b tt rec_surr z1 env_trg_pa env_trg_npa env_trg_surr name
    
    ctx=ecog(j,:);
    %%% ecog_shuf=suffle_data(ctx);   optional surrogate
    Ecogfiltered=filtfilt(b,a,ctx);
    
    data_ones=find(data(j,:)==1);
    hp=wrapToPi(angle(hilbert(Ecogfiltered)));
    ang=hp(data_ones);
    if (circ_rtest(ang))<0.05 
        
        %%% clear Ecogfiltered;Ecogfiltered=filtfilt(b,a,ecog_shuf);  optional surrogate
        env=abs(hilbert(Ecogfiltered));
        
        data_g=smoothdata(data(j,:),'gaussian',25); % convolve a gaussian to spikes
        %%%option: N = 100; alpha =15 ; w = gausswin(N,alpha); data_g=conv(data(j,:),w,'same');
                
        onset1=bursts(env);
        onset1=horzcat(onset1{:}); % onset of short and long bursts
        onset=bursts_aligned(env,Ecogfiltered); 
        onset=horzcat(onset{:}); % onset of short and long bursts aligned to closer peak of beta oscillations
        
        epoch=200;
        for n=1:(length(env)/100)
            idx_sur=randi([epoch+1,(length(env)-epoch)],1,1);
            pre_sur(n,:)= data_g(idx_sur-epoch:idx_sur+epoch);
        end
        
        
        for jj=1:size(onset,2)
            if onset(jj)>200 && onset(jj)+200<length(data_g) && onset1(jj)>200 && onset1(jj)+200<length(data_g)
                output_pa(jj,:)=data_g(onset(jj)-200:onset(jj)+200);
                output_npa(jj,:)= data_g(onset1(jj)-200:onset1(jj)+200);
            end
        end
        
        if max(zscore(sum(output_pa,1)))>1.96
            tt=tt+1;
            count_b(1,tt)=length(onset);
            rec_pa(tt,:)=zscore(sum(output_pa,1));
            rec_npa(tt,:)=zscore(sum(output_npa,1));
            rec_surr(tt,:)=zscore(sum(pre_sur,1));
            
            %envelope of filtered triggred_avg for statistical comparison
            env_trg_pa(tt,:)=abs(hilbert(filtfilt(b,a,rec_pa(tt,:))));
            env_trg_npa(tt,:)=abs(hilbert(filtfilt(b,a,rec_npa(tt,:))));
            env_trg_surr(tt,:)=abs(hilbert(filtfilt(b,a,rec_surr(tt,:))));
            
        end
    end
end


region_pl=(mean(rec_pa,1));
region_spl=(std(rec_pa)./sqrt(size(rec_pa,1)));
region_npl=(mean(rec_npa,1));
region_snpl=(std(rec_npa)./sqrt(size(rec_npa,1)));
region_surr=(mean(rec_surr,1));
region_ssurr=(std(rec_surr)./sqrt(size(rec_surr,1)));


%%%plot envelopes 
% figure()
% subplot(1,4,1)
% plot(env_trg_pa')
% subplot(1,4,2)
% plot(env_trg_npa')
% subplot(1,4,3)
% plot(env_trg_surr')
% subplot(1,4,4)
% plot(mean(env_trg_pa))
% hold on
% plot(mean(env_trg_npa))
% plot(mean(env_trg_surr))

time2=1:401;
color_s=[0.5 0.5 0.5];
load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','squash','blood','sky','aegean');
if name=='bz'
    color_b=blood;
else
    color_b=aegean;
end


%%% plot of phase aligned trig-avg plus stat significance with non-phase
%%% aligned trig_avg and surrogates
fig1=figure()
y2=region_pl; y1=y2+region_spl; y3=y2-region_spl;
y5=region_npl; y4=y5+region_snpl; y6=y5-region_snpl;
p1=plot(time2, y2, 'LineWidth',1.5,'Color',color_b)
patch([time2 fliplr(time2)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.1],'EdgeColor','none')
patch([time2 fliplr(time2)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.1],'EdgeColor','none')
hold on
for mm=1:2
    if mm==1
        clear A; A=env_trg_pa;
        clear B; B=env_trg_npa;
        color_c=[0.5 0.5 0.5];
        place=[2.9 2.9 3.0 3.0];
    else
        clear A; A=env_trg_pa;
        clear B; B=env_trg_surr;
        color_c=[0 0 0];
        place=[3.3 3.3 3.4 3.4];
    end
    time=1:401;
    st=NaN(1,401);
    hayriye_c; st(1,:)=stats.prob; st2(1,:)=stats.posclusterslabelmat;
    beg=find(st(1,:)<0.01 & st2(1,:)~=0);
    if ~isempty(beg)
        sig_rise_all=[beg(1) beg(end)];
    else
        sig_rise_all=[];
    end
    if ~isempty(sig_rise_all)
        sig{1,mm}=sig_rise_all;
%         patch([sig_rise_all(1) sig_rise_all(2) sig_rise_all(2) sig_rise_all(1)],place,color_c,'EdgeColor','none')
    else
        sig{1,mm}=[];
    end
    
end
xline(200,'--',{'burst onset'},'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','Color',[0.5 0.5 0.5],'LineWidth',2)
xlim ([0 400])
ylim ([-3 3])
yticks([-3:1:3])
xticks([0:100:400])
xticklabels ({'-200','-100','0','100','200'})
box ('off')
xlabel ('Time (ms)')
ylabel({'Firing-rate','(z-score)'})
set(gca,'FontSize',12)

fig1.Units = 'centimeters';
fig1.OuterPosition= [10, 10, 8, 10];
fig1.Color='w';




load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','squash','blood','sky','aegean','grey','black');
if name=='bz'
    col={blood,grey,black};
else
    col={aegean,grey,black};
end
cond={'env_trg_pa','env_trg_npa','env_trg_surr'};
label={'phase-aligned','non phase-aligned','random ref'};
fig1=figure
for ii=1:size(cond,2)
    data=eval(cond{1,ii});
    color_b=col{1,ii};
    y2=mean(data);
    y1=y2+(std(data)./sqrt(size(data,1)));
    y3=y2-(std(data)./sqrt(size(data,1)));
    p1=plot(time2, y2, 'LineWidth',1.5,'Color',color_b)
    patch([time2 fliplr(time2)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.1],'EdgeColor','none')
    patch([time2 fliplr(time2)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.1],'EdgeColor','none')
hold on
end

place={[2.5 2.5 2.6 2.6] [2.8 2.8  2.9 2.9]};
for mm=1:2
if ~isempty(sig{1,mm})
    patch([sig{1,mm}(1) sig{1,mm}(2) sig{1,mm}(2) sig{1,mm}(1)],place{1,mm},col{1,mm+1},'EdgeColor','none')
end
end

xline(200,'--',{'burst onset'},'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','Color',[0.5 0.5 0.5],'LineWidth',2)
xlim ([0 400])
ylim ([-3 3])
yticks([-3:1:3])
xticks([0:100:400])
xticklabels ({'-200','-100','0','100','200'})
box ('off')
xlabel ('Time (ms)')
ylabel({'envelope \beta filtered','firing-rate (zscore)'})
set(gca,'FontSize',12)

fig1.Units = 'centimeters';
fig1.OuterPosition= [10, 10, 8, 10];
fig1.Color='w';


%%% suplots: phase_aligned trig-acg, non-phase aligned avg and surrogates
fig=figure;
subplot(1,3,1)
y2=region_pl; y1=y2+region_spl; y3=y2-region_spl;
p1=plot(time2, y2, 'LineWidth',1.5,'Color',color_b)
patch([time2 fliplr(time2)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.1],'EdgeColor','none')
patch([time2 fliplr(time2)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.1],'EdgeColor','none')
xline(200,'--',{'burst onset'},'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','Color',[0.5 0.5 0.5],'LineWidth',2)
xlim ([0 400])
ylim ([-3 3])
yticks([-3:1:3])
xticks([0:100:400])
xticklabels ({'-200','-100','0','100','200'})
legend([p1],{'phase-aligned'},'FontSize',12,'box','off','Location','northeast' )
box ('off')
xlabel ('Time (ms)')
ylabel('Firing-rate(z-score)')
set(gca,'FontSize',12)


subplot(1,3,2)
y5=region_npl; y4=y5+region_snpl; y6=y5-region_snpl; color_s=[0.5 0.5 0.5];
p2=plot(time2, y5, 'LineWidth',1.5,'Color',color_s)
patch([time2 fliplr(time2)], [y4 fliplr(y5)],[color_s],'FaceAlpha',[0.1],'EdgeColor','none')
patch([time2 fliplr(time2)], [y5 fliplr(y6)],[color_s],'FaceAlpha',[0.1],'EdgeColor','none')
xline(200,'--',{'burst onset'},'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','Color',color_s,'LineWidth',2)
xlim ([0 400])
ylim ([-5 5])
xticks([0:100:400])
xticklabels ({'-200','-100','0','100','200'})
legend([p2],{'non-phase aligned'},'FontSize',12,'box','off','Location','northeast' )
box ('off')
xlabel ('Time (ms)')
ylabel('Firing-rate(z-score)')
set(gca,'FontSize',12)

subplot(1,3,3)
y5=region_surr; y4=y5+region_ssurr; y6=y5-region_ssurr; color_s=[0 0 0];
p2=plot(time2, y5, 'LineWidth',1.5,'Color',color_s)
patch([time2 fliplr(time2)], [y4 fliplr(y5)],[color_s],'FaceAlpha',[0.1],'EdgeColor','none')
patch([time2 fliplr(time2)], [y5 fliplr(y6)],[color_s],'FaceAlpha',[0.1],'EdgeColor','none')
xline(200,'--',{'burst onset'},'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','Color',color_s,'LineWidth',2)
xlim ([0 400])
ylim ([-5 5])
xticks([0:100:400])
xticklabels ({'-200','-100','0','100','200'})
legend([p2],{'surrogate'},'FontSize',12,'box','off','Location','northeast' )
box ('off')
xlabel ('Time (ms)')
ylabel('Firing-rate(z-score)')
set(gca,'FontSize',12)

fig.Units = 'centimeters';
fig.OuterPosition= [10, 10, 30, 10];
fig.Color='w';

end