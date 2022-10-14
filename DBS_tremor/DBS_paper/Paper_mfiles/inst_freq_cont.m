
function[dev_freq]=inst_freq_cont(s,start,ending,tp_s,tp_e,match_ax,pt)

%%%% you need 1 second (1000 ponits segments)

peaki=[NaN 4 NaN 3.5];

dum2=s.freq{1,3};
freqi=dum2(match_ax(2,pt,1),:);
ch_fr=abs(freqi-peaki(pt));
dev_freq =[mean(ch_fr) std(ch_fr)];

[ref]=pt_post(dum2, start, ending,tp_s,tp_e); %%% uncomment plot to see breaks in data


for t=1:size(ref,1)
    r=0;
    for j = 2:2:length(ref{t,1})
        r=r+1;
        if ~isempty(length(ref{t,1}))
            tempo=floor((numel(ref{t,1}(j-1):ref{t,1}(j))/1000))*1000;
            period=ref{t,1}(j)-tempo:ref{t,1}(j);
            bins=1:1000:length(period);
            for g=1:length(bins)
                if (g+1)<length(bins)
                        y(1,g) = mean(dum2(match_ax(2,pt,1),period(bins(g):bins(g+1))));
                end
            end
        end

        
        n=abs(y-peaki(pt));
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
        total(r,:)=[numel(find(p==1))./length(p) numel(find(p==2))./length(p) numel(find(p==3))./length(p) numel(find(p==4))./length(p)];
        end
        
        clearvars -except s samplerate ref color_b1 match_ax total t j dum dum2 r pt peaki
    end




    f1=figure(t+1)
    for u=1:3
        subplot(1,3,u)
        pie(total(u,:))
    end
    box('off')
    f1.Units = 'centimeters';
    f1.OuterPosition= [10, 10, 35, 8];
    set(f1,'color','w')

    clear total
end

    legend('{<0.5HZ}','{0.5-1Hz}','{1-2Hz}','{>2Hz}')
    legend('boxoff')
    title(sprintf('patient %d',(pt)))

end


