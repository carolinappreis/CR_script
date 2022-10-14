
    close all, clear all
    load('/Users/Carolina/Desktop/oxford/data_code_thesis/DBS/paper_analyses/paper_matfiles/instantaneous_freq_rs.mat')    %%% instantaneous frequency in posture and spiral , 3 axis plus peak freq and the main axis (from DBS_newmaster)
co=2;

figure(1)
title('Posture RSS')
figure(2)
title('Spiral RSS')

peaki=[4.5 4 4 3.5];

for spiral1=1:2
    
    for iii=1:4
    clearvars -except spiral1 iii cf co
        spiral=spiral1-1; [match_ax]=link_ax(spiral);

        paxis_frq=cf.ifrq{iii,spiral1}(match_ax(co,iii,1),cf.h_up{iii,spiral1});

        c_freq=peaki(iii);
        y=[];

        bins=1:1000:length( paxis_frq);
        for g=1:length(bins)
            if (g+1)<length(bins)
                ax=1;
                %                 for ax=1:size(dum,1)
                y(ax,g)= mean(paxis_frq(1,bins(g):bins(g+1)));
                %                 end
            end
        end

        n=abs(y-c_freq);
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

        total=[numel(find(p==1))./length(p) numel(find(p==2))./length(p) numel(find(p==3))./length(p) numel(find(p==4))./length(p)];
       
        p1=figure(spiral1);
        subplot(1,4,iii)
        pie(total)
        set(p1,'color','w');
        title(sprintf('patient %d',(iii)))
        

    end
    legend('{<0.5HZ}','{0.5-1Hz}','{1-2Hz}','{>2Hz}')
    legend('boxoff')
end