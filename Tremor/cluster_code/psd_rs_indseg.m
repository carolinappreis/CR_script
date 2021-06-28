function [out]=psd_rs_indseg(out)
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/pxx_ns_indivseg.mat')
pNSSA(1,:)=pxx_nsseg;

samplerate=1000;
for iii=1:size(out.start_c,1)
    m_ax=1;
    ARC=squeeze(out.change_c{iii,2}{m_ax,1});
    sup1=[];
    amp1=[];
    
    for y1=1:size(ARC,2)
        for y2=1:size(out.z_seg{iii,y1},1)
            dum=(out.z_seg{iii,y1}(y2,:));
            [pp,F]=pwelch(dum,samplerate,[],samplerate,samplerate);
            frange=find(F==3):find(F==8);
            pxx_maxf=max(pp(frange,1));
            
            
            if ARC(y2,y1)>0
                amp1=[amp1   pxx_maxf];
            else
                sup1=[sup1   pxx_maxf];
            end
            clear  pxx_maxf dum pp F frange
        end
    end
    
    
    pNSSA(2,iii)=nanmean(sup1);
    pNSSA(3,iii)=nanmean(amp1);
    
    clear amp1 sup1   
end



pNSSA(:,3)=NaN;


ttest(pNSSA(1,:), pNSSA(2,:))
ttest(pNSSA(1,:), pNSSA(3,:))
ttest(pNSSA(2,:), pNSSA(3,:))

load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','aegean','stone')
color_b1={ stone; aegean; blushred};


bar(1:3,[nanmedian(pNSSA(1,:)) nanmedian(pNSSA(2,:)) nanmedian(pNSSA(3,:))],'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.3,'EdgeColor','none')
hold on
for i=1:size(pNSSA,2)
    cl= rand(1,3);
    plot(pNSSA(1:3,i),'Color',cl,'LineWidth',1)
    plot([pNSSA(1,i); pNSSA(2,i); pNSSA(3,i)],'.','MarkerSize',15,'Color',cl)
    xlim([0 4])
    xticklabels({'NS','SUP','AMP'})
    box('off')
    ylabel({'Power at tremor' ; 'frequency peak (u.a)'})
    set(gca,'FontSize',12)
end



% ttest(pNSSA(2,:), pxx_ns_sup)
% ttest(pNSSA(3,:), pxx_ns_amp)
% 
% bar(1:2,[[median(pNSSA(2,:)) ; median(pNSSA(3,:))] [ median(pxx_ns_sup); median(pxx_ns_amp)]]);
% 
% data={[pNSSA(2,:)' pNSSA(3,:)'] , [pxx_ns_sup' pxx_ns_amp']};
% boxplotGroup(data, 'PrimaryLabels', {'stim' 'nostim'},'SecondaryLabels',{'SUP', 'AMP'}, 'GroupLabelType', 'Vertical')

end

