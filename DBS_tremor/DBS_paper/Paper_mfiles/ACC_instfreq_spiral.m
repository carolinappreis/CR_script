cd('/Users/Carolina/Documents/GitHub/CR_script/DBS_tremor/DBS_paper')

clear all
ad=struct; ad.metrics=cell(1,1);ad.metricnames={'subj';'cond';'trials';'mean env';'peak2peak';'variance'};

load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/ACC_zscore_acr_cond.mat'))

%%% normalised to NS condition: load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/ACC_spiral_NSnorm.mat'))


% (smooth((1000/(2*pi))*diff(unwrap(s.phase{iii,co}(i,:))),500))';

peaki=[4 4 3.5];

for iii=1:3
    clearvars -except ad iii signal_all map peaki total total2 dev_freq
    for co=3  %%% co={'NS';'HF';'C'};
        pha=angle(hilbert(signal_all{iii}));
        freqi=(smooth((1000/(2*pi))*diff(unwrap(pha)),500))';
        dev_freq(iii,:)=[mean(abs(freqi-peaki(iii))) std(abs(freqi-peaki(iii)))];

        r=0;


        for trial= 1:size(map{iii,co},1)
            r=r+1;
            period=map{iii,co}(trial,1):map{iii,co}(trial,end);
            bins=1:1000:length(period);
            for g=1:length(bins)
                if (g+1)<length(bins)
                    y(1,g) = mean(freqi(1,period(bins(g):bins(g+1))));
                end
            end
            clear period bins



            n=abs(y-peaki(iii));
            if max(n)<4
                p=[];
                for ii=1:length(n)
                    if n(ii)<0.5
                        p(1,ii)=1;
                    elseif n(ii)>=0.5 && n(ii)<1
                        p(1,ii)=2;
                    elseif n(ii)>=1 && n(ii)<2
                        p(1,ii)=3;
                    else
                        p(1,ii)=4;
                    end
                end
            else
                error('outlier')
            end


            for ii=1:size(p,1)
                total{iii,1}(r,:)=[numel(find(p==1))./length(p) numel(find(p==2))./length(p) numel(find(p==3))./length(p) numel(find(p==4))./length(p)];
            end

        end
    end
    total2(iii,:)=mean(total{iii,1});

    f1=figure(1)
    subplot(1,3,iii)
    pie(total2(iii,:))
    box('off')
    f1.Units = 'centimeters';
    f1.OuterPosition= [10, 10, 35, 8];
    set(f1,'color','w')
   title(sprintf('patient %d',(iii+1)))



end

    legend('{<0.5HZ}','{0.5-1Hz}','{1-2Hz}','{>2Hz}')
    legend('boxoff')



   t_sup=[0.78 ; 1.66 ; 0.315];

  

[fitresult] = fit(t_sup, dev_freq(:,1),'poly1' );
h=plot(fitresult,t_sup, dev_freq(:,1));
legend('off')
hold on
plot(t_sup, dev_freq(:,1),'k.','MarkerSize',10,'HandleVisibility','off');
[r,p]=corrcoef(dev_freq(:,1),t_sup)
xlim([0 2])
xlabel('tremor supression (NS-PSS)')
ylabel('absolute deviation of freq (stim freq - inst freq)')
box('off')