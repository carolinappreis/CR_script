function [h]=twitch(iii,d,samplerateold,samplerate,co)

if co==2
    data=d.data_raw;
    
    %%% determine stimulation time points
    index=[];
    for i=2:size(data,2)-1
        if data(2,i-1)<2.5 && data(2,i)>2.5
            index=[index i];
        end
    end
    
    main=[1 1 3 1 3 3 3 3 1 1];
    dum=[3 5 6];
    
    zdata=zscore(d.data_raw(dum(main(iii)),:));
    
    f=floor(((250/2)./samplerate)*samplerateold);
    seg_raw=NaN(length(index),f*2);
    seg_pha=NaN(length(index),f*2);
    
    for i=1:length(index)
        if (~isnan(index(i)))
            seg_raw(i,1:(f*2))=zdata(1,(index(i)-(f-1)):(index(i)+f));
        else
            seg_raw(i,1:(f*2))=NaN;
        end
    end
    
    figure(1)
    data1=seg_raw;
    subplot(2,5,iii)
    y2=nanmean(data1,1);
    sem = nanstd(data1,1)./ sqrt(size(data1,2));
    y1=y2+sem;
    y3=y2-sem;
    time=1:length(y2);
    plot(y2,'k','LineWidth',1);
    hold on
    patch([time fliplr(time)], [y1 fliplr(y2)],'k','FaceAlpha',[0.15],'EdgeColor','none','HandleVisibility','off')
    patch([time fliplr(time)], [y2 fliplr(y3)],'k','FaceAlpha',[0.15],'EdgeColor','none','HandleVisibility','off')
    hold on
    xline(size(data1,2)/2,'--')
    box('off')
    xlim([0 size(data1,2)])
%     ylim ([-1 1])
    xticks([size(data1,2)/2])
    xticklabels(['stim'])
    yticks([])
    title(sprintf('patient %d',(iii)))
    
%     
%     if iii==9 || iii==10  %%% the only patients with EMG on the thenar eminence
%         
%         zdata1=zscore(d.data_raw(9,:)); %emg 3 for patient 9 and 10
%         for i=1:length(index)
%             if (~isnan(index(i)))
%                 seg_raw1(i,1:(f*2))=zdata1(1,(index(i)-(f-1)):(index(i)+f));
%             else
%                 seg_raw1(i,1:(f*2))=NaN;
%             end
%         end
%         
%         figure(2)
%         subplot(2,1,iii-8)
%         data1=seg_raw1;
%         y2=nanmean(data1,1);
%         sem = nanstd(data1,1)./ sqrt(size(data1,2));
%         y1=y2+sem;
%         y3=y2-sem;
%         time=1:length(y2);
%         plot(y2,'k','LineWidth',1);
%         hold on
%         patch([time fliplr(time)], [y1 fliplr(y2)],'k','FaceAlpha',[0.15],'EdgeColor','none','HandleVisibility','off')
%         patch([time fliplr(time)], [y2 fliplr(y3)],'k','FaceAlpha',[0.15],'EdgeColor','none','HandleVisibility','off')
%         hold on
%         xline(size(data1,2)/2,'--')
%         box('off')
%         xlim([0 size(data1,2)])
%         ylim ([-1 1])
%         xticks([size(data1,2)/2])
%         xticklabels(['stim'])
%         yticks([])
%         title(sprintf('patient %d',(iii)))
%         
%     end
    h=1;
end
h=2;

end





