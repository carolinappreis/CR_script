% clear; close
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/cluster_out_mc.mat','out');
% m_ax=1;% change if the main axis is not always 1 - replace by array of main axes
% load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','aegean','stone','squash','sapphire','azure');
% color_b1=[blushred; aegean; stone];


for iii=10
%     1:size(out.start_c,1)
    %%% for type=1:2
    close all
    type=1; feature={'out.change_c';'out.fchange'};
    dr=eval([(num2str(feature{type,1})) '{' num2str(iii) ',2}{m_ax,1}']);
    
    data=nanmedian(dr);
    nostim=eval(squeeze([num2str(feature{type,1}) '{' num2str(iii) ',1}(m_ax,:)']));
    n=[];
    n=[n ; data(find(data > prctile(nostim, 99.7917) | data< prctile(nostim, 0.2083)))];
    n_phases=numel(find(data > prctile(nostim, 99.7917) | data< prctile(nostim, 0.2083)));
    if ~isempty(find(data > prctile(nostim, 99.7917) | data< prctile(nostim, 0.2083)))
        phase=find(data > prctile(nostim, 99.7917) | data< prctile(nostim, 0.2083));
    else
        phase=NaN;
    end
    
    
    
    
    f1=figure;
    for ii=1:size(gg,1)
        subplot(size(gg,1),1,ii)
        y1=smooth(gg(ii,:));
        x1=1:length(y1);
        plot((y1),'Color','k','LineWidth',2)
        hold on
        mdp1=fitlm(x1,y1,'poly2');
        plot(mdp1.Fitted,'r','LineWidth',1.5)
        n=mdp1.Rsquared.Adjusted;
        
        filename=['adjr^2=',num2str(round(n,2))];
        t = text(0.5,0.5,filename);
        s = t.FontSize;
        t.FontSize = 12;
        
        xline(ini,'--','Color',[0.5 0.5 0.5],'LineWidth',2)
        xlim([0 80000])
        xticks([0:10000:80000])
        xticklabels({'-20','-10','0','10','20','30','40','50','60'})
        ylim([0 14])
        yticks(0:2:14)
        ylabel('tremor severity (m/s^2)')
        xlabel ('time (s)')
        

        box('off')
        set(gca,'FontSize',12)
        title(['PLS pt',num2str(numb),' trial',num2str(ii)])
        
    end
    
    f1.Units = 'centimeters';
    f1.OuterPosition= [10, 10, 14, 16];
    set(f1,'color','w');
    box('off')
    
    
    
    
    
    if (~isnan(phase))
        for p=1:length(phase)
            dat=out.zenv_seg{iii,phase(p)};
            f1=figure(p)
            for gg=1:size(dat,1)
                subplot(2,5,gg)
                x=1:size(dat,2);
                y=dat(gg,:);
                plot(x,y,'k')
                hold on
                mdp1=fitlm(x,y,'poly1');
                plot(mdp1.Fitted,'r','LineWidth',1.5)
                n=round(mdp1.Rsquared.Adjusted,2);
                title(sprintf('r= %d',(n)))
            end
            set(f1,'color','w');
            f1.OuterPosition= [1,339,1440,539];
        end
        
    end
    
    
end

yy(numb,:)=smooth(median(gg,1));
%     y=yy(numb,1:35000);
y=yy(numb,1:40000);
x=tt(1:length(y));
initial_params=[];
[param]=sigm_fit(x,y,initial_params)        % automatic initial_params
clear x y