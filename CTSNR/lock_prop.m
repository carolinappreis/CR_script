  function [fig]=lock_prop(data,ecog,name)
freq={[1 7];[8 12];[15 35];[36 80];[81 100]};
srn=1000;
for t=1:size(freq,1)
    for r=1:size(data,1)
            dati=data(r,:);
            sk_ctc=[];
            for i =1:srn:(length(dati)-srn);
                sk_ctc=[sk_ctc numel(find(dati(i:i+srn)==1))];
            end
            data_ones=find(dati==1);
            [b,a]=butter(2,[(freq{t,1}(1))/(0.5*srn) (freq{t,1}(2))/(0.5*srn)],'bandpass');
            Ecogfiltered=filtfilt(b,a,ecog(r,:));
            ang=angle(hilbert(Ecogfiltered(data_ones)))';
            vec_an(r,:)=circ_r(ang);
            vec_pref(r,:)=circ_mean(ang);
            spikerate(r,:)=mean(sk_ctc);
            clear ang Ecogfiltered b a data_ones dati;
    end
    vec(t,:)=mean(vec_an);
    pref(t,:)=circ_mean(vec_pref);
    pref_std(t,:)=circ_std(vec_pref);
    clear vec_an vec_pref dati 
end




fig=figure()
subplot(1,3,1)
histogram(spikerate)
ylabel('Unit count');
xlabel('Frequency (Hz)');
title('Firing rate')
box('off')

subplot(1,3,2)
bar(vec')
title ('Phase locking of units to ECoG')
xticklabels({'1-7','8-12','15-35','36-80','81-100'})
ylabel('Vector length');
xlabel('Frequecy (Hz)');
box('off')

subplot(1,3,3)
i=3;
polarplot([0 pref(i)], [0, vec(i)],'linewidth',2)
hold on
polarplot([0 pref(i)+pref_std(i)], [0, vec(i)],'linewidth',2,'Color',[0.5 0.5 0.5])
polarplot([0 pref(i)-pref_std(i)], [0, vec(i)],'linewidth',2,'Color',[0.5 0.5 0.5])
title ('Preferred phase of units to \betaECoG')
fig.Color='w';
box('off')
fig.Units = 'centimeters';
fig.OuterPosition= [10, 10, 36, 10];
fig.Color='w';



  end