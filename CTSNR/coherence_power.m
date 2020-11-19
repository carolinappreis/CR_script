clear all
close all
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/final_mats')
load('SNR_bua.mat');

c=0;
p=0;
samprate=1000;
for pr=1:size(bua,1)
    c=c+1;
    for ct=2:size(bua{pr,1},1)
        [Pxx_ind,F_i]=mscohere(bua{pr,1}(1,:),bua{pr,1}(ct,:),samprate,[],samprate,samprate);
        frange=find(F_i==15):find(F_i==35);
        if (sum(Pxx_ind(frange))/sum(Pxx_ind(1:end)))>0.1
            p=p+1;
            coher(p,:)=Pxx_ind;clear Pxx_i
            [Pxx_ctx,F_ind]=pwelch(bua{pr,1}(1,:),samprate,[],samprate,samprate);
            power_ctx(c,:)=Pxx_ctx; clear Pxx_ctx
            [Pxx_probe,F_ind]=pwelch(bua{pr,1}(ct,:),samprate,[],samprate,samprate);
            power_probe(p,:)=Pxx_probe; clear Pxx_probe
        end
    end
end

load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','squash','blood','sky','aegean');
if name=='bz'
    color_b={squash blood};
else
    color_b={sky aegean};
end

plt={'power_ctx';'power_probe';'coher'};
label={'Log Cortical Power' ; 'Log Subcortical Power';'Mean Squared Coherence';};




fig=figure;
 time=F_i';
for ii=1:size(plt,1)
    subplot(1,size(plt,1),ii)
    data=eval(plt{ii,1});
    if ii==3
        y2=mean(data);
        y1=y2+std(data)./sqrt(size(data,1));
        y3=y2-std(data)./sqrt(size(data,1));
       
        
        box('off')
    else
        y2=log(mean(data));
        y1=log(mean(data)+std(data)./sqrt(size(data,1)));
        y3=log(mean(data)-std(data)./sqrt(size(data,1)));
    end
        plot(time, y2,'LineStyle','-', 'LineWidth',1.5,'Color',[color_b{1,2}])
        hold on
        patch([time fliplr(time)], [y1 fliplr(y2)],[color_b{1,2}],'FaceAlpha',[0.2],'EdgeColor','none')
        patch([time fliplr(time)], [y2 fliplr(y3)],[color_b{1,2}],'FaceAlpha',[0.2],'EdgeColor','none')
        xlim ([0 80])
        xticks([0:10:80])
        box('off')
        fig.Units = 'centimeters';
        fig.OuterPosition= [10, 10, 30, 10];
        fig.Color='w';
        set(gca,'FontSize',12)
        ylabel(sprintf(label{ii,1}));
        xlabel('Frequency (Hz)');
end