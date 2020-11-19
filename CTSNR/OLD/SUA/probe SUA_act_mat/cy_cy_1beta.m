clear all
%  cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\probe SUA_act_mat')
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/probe SUA_act_mat')
load('data_SUA_SNR.mat')

srn=1000;
[b,a]=butter(2,[15/(0.5*srn) 35/(0.5*srn)],'bandpass');
data=cell2mat(data_region);
ecog=[];
for i=1:size(data_region,1)
    ecog = [ecog ; repmat(Ecog_region(i,:),size(data_region{i,1},1),1)];
end


tt=0;
for j=1:size(data,1)
    
    clearvars -except j b a ecog data srn tt  check vec_lg pref_pha
    
    ctx=ecog(j,:);
    ecog_shuf=suffle_data(ctx);
    
    Ecogfiltered=filtfilt(b,a,ctx);
    %      Ecogfiltered=filtfilt(b,a,ecog_shuf);
    
    env=abs(hilbert(Ecogfiltered));
    
    data_ones=find(data(j,:)==1);
    hp=wrapToPi(angle(hilbert(Ecogfiltered)));
    ang=hp(data_ones);
    
    if (circ_rtest(ang))<0.05
        
        cy_bursts=cycles_10(env,Ecogfiltered);
        block=cell2mat(cy_bursts);

        for d1=1:size(block,1)
            for d2=1:size(block,2)
                if d2+1<length(block(d1,:))
                    epoch=block(d1,d2):block(d1,d2+1);
                    l=find(data(j,epoch)==1);
                    if ~isempty(epoch(l))
                        pha_b(d1,d2)=hp(epoch(l(1))); %%% picking just the first spike in a cycle l(1)
                    else
                        pha_b(d1,d2)=NaN;
                    end
                end
            end
        end
        
       
        for x=1:size(pha_b,2)
            bu=pha_b(:,x); bu=bu(~isnan(bu));
            check1(1,x)=length(bu);
%             if ~isempty(bu) && length(bu)>=20
%                 bu1=bu(randperm(20));
%             else
%                 bu1=NaN;
%             end
            vl(1,x)=circ_r(bu);
            pp(1,x)=circ_mean(bu); clear bu bu1
        end
        
        if min(check1)>size(block,1)/5
        tt=tt+1;
        vec_lg(tt,:)=vl;
        pref_pha(tt,:)=pp;
        check(tt,:)=check1;
        end
    end
end


%%% cycle eleven is burst onset

fig=figure()
cy_plots=7:15;
for i=1:length(cy_plots)
    subplot(1,length(cy_plots),i,polaraxes)
    for un=1:size(vec_lg,1)
        polarplot([0 pref_pha(un,cy_plots(i))], [0, vec_lg(un,cy_plots(i))],'linewidth',2)
        rlim([0 0.5])
        hold on
    end
end

fig.Units = 'centimeters';
fig.OuterPosition= [10, 10, 40, 10];
fig.Color='w';



vec_m=nanmean(vec_lg,1);
ang_m=circ_mean(pref_pha);

fig=figure()
cy_plots=6:16; % first 5
for i=1:length(cy_plots)
    if ~isempty(intersect(cy_plots(i), [1:10]))
        cl=[0.5 0.5 0.5];
    elseif cy_plots(i)==11
        cl=[0.5 0 0];
    else
        cl=[0 0 0.5];
    end
    if i==1
        polarplot([0 ang_m(1,cy_plots(i))], [0, vec_m(1,cy_plots(i))],'Color',cl,'linewidth',1,'MarkerIndices',[3],'Marker','d')
    else
        polarplot([ang_m(1,cy_plots(i-1)) ang_m(1,cy_plots(i))], [vec_m(1,cy_plots(i-1)), vec_m(1,cy_plots(i))],'Color',cl,'linewidth',1,'MarkerIndices',[2],'Marker','d')
    end
    hold on
end

color_b=[0 0 0.5];
fig=figure()
% plot(nanmean(vec_lg,1),'-d','LineWidth',1.5,'Color',color_b)
plot(mean(vec_lg),'-d','LineWidth',1.5,'Color',color_b)
hold on
xline(11,'--',{'burst onset'},'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','Color',[0.5 0.5 0.5],'LineWidth',2)
xlim([0 22])
xlabel('Number of {\beta} cycles')
xticks([1:2:21])
xticklabels ({'-10','-8','-6','-4','-2','0','2','4','6','8','10'})
ylabel('Vector length')
box('off')
set(gca,'FontSize',12)
fig.Units = 'centimeters';
fig.OuterPosition= [10, 10, 10, 10];
fig.Color='w';
