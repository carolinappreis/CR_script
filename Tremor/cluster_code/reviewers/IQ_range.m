% close all; clearvars -except out
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/bini.mat');
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/cluster_out_mc.mat','out');

load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','stone');
cl=stone;
m_ax=1;
for iii=1:10
    th=0.5;
    f1=figure(14)
    subplot(2,5,iii)
    y=amp_n_bins(iii,:);
    x=bins(1:end-2);
    r=find(y<th);
    dum=find(diff(r)~=1);
    
    if iii==9
        a=2.25;
        b=5;
    else
        
        if ~isempty(dum) && numel(r)==1 && r>numel(x)/2
            [val,idx]=min(abs(y(1:round((length(y)/2)))-th));
            b=x(r);
            a=(x(idx)-1);
            
        elseif ~isempty(dum) && numel(r)==1 && r<numel(x)/2
            [val,idx]=min(abs(y((round(length(y)/2)):end))-th);
            b=x(r);
            a=(x(idx+round(length(y)/2)));
            
        elseif  ~isempty(dum) && numel(r)>1
            b=x(r(dum+1));
            a=x(r(dum));
            if iii==5
%                 b=x(r(dum+1)+2);
                a=5;b=7;
                
            end
        else
            q=find(isnan(y));
            w= find(diff(q)~=1) ;
            if ~isempty(w)
                b=x(q(w+1))
                a=x(q(w));
            else
                if q(1)>length(x)/2
                    b=x(q(1)); a=x(1)
                else
                    b=x(end); a=q(end);
                end
            end
            
        end
    end
    
    range(iii,1)=b-a;
    r_numb(iii,1:2)=[a b];
    bar(x,y,'FaceColor',cl,'EdgeColor',cl)
    hold on
    yline(th)
    xline(a,'r')
    xline(b,'r')
    xticks([1:1:14])
    %     xticklabels({'2','3','4','5','6','7','8'});
    ylabel({'Absolute amplitude (\muV^2)'})
    xlabel('Frequency (Hz)')
    %     set(gca,'FontSize',12)
    box('off')
    legend('off')
    
    clear y full
    f1.Units = 'centimeters';
    f1.OuterPosition= [10, 10, 25, 10];
    set(f1,'color','w');
    
    d1=nanmedian(squeeze(out.change_c{iii,2}{1,m_ax}));
    arc.value(iii,1,1)=nanmean(abs(d1))*100;
    arc.value(iii,2,1)=(max(d1)-min(d1));
    if min(d1)<0
        arc.value(iii,3,1)=-(min(d1)*100);
    else
        arc.value(iii,3,1)=NaN;
    end
    %     arc.value(iii,4,1)=max(d1)*100;
    arc.label={'mean effect';'effect range';'supressive effect'};
    clear d1
end
clear x
% for ii=1
%     
%     f1=figure(15)
%     
%     
%     dum=find(~isnan((arc.value(:,ii))));
%     x=range(dum);
%     
%     y=arc.value(dum,ii);
%     
%     [fitresult] = fit( x, y, 'poly1' );
%     h=plot(fitresult,x,y);
%     
%     hold on
%     plot(x,y,'k.','MarkerSize',10,'HandleVisibility','off');
%     legend('off')
%     [r(1,ii),p(1,ii)]=corrcoef(x,y);
%     
%     %     legend(h,['r=',num2str(r)],'box','off')
%     box('off')
%     
%     xlabel('IQ range')
%     ylabel(sprintf( arc.label{ii,1}))
%     
%     
%     set(gca,'FontSize',12)
%     f1.Units ='centimeters';
%     %     f1.OuterPosition= [10, 10, 8, 10];
%     set(f1,'color','w');
%     
%     %     y2=plot(x(sig),y(sig),'ko','MarkerSize',8);
%     
% end



%%%% run f=if th=0.5

r_numb(2,2)=5.25;
r_numb(4,2)=6.25;
r_numb(6,1)=4.25;
r_numb(7,2)=4.5;
r_numb(9,1)=3.5;
r_numb(9,2)=4.75;
r_numb(10,2)=6.25;

for iii=1:10
    
     th=0.5;
    f1=figure(14)
    subplot(2,5,iii)
range(iii,1)=r_numb(iii,2)-r_numb(iii,1);
y=amp_n_bins(iii,:);
    x=bins(1:end-2);
    bar(x,y,'FaceColor',cl,'EdgeColor',cl)
    hold on
    yline(th)
    xline(r_numb(iii,1),'r')
    xline(r_numb(iii,2),'r')
    xticks([1:1:14])
    %     xticklabels({'2','3','4','5','6','7','8'});
    ylabel({'Absolute amplitude (\muV^2)'})
    xlabel('Frequency (Hz)')
    %     set(gca,'FontSize',12)
    box('off')
    legend('off')
    
    clear y full
    f1.Units = 'centimeters';
    f1.OuterPosition= [10, 10, 25, 10];
    set(f1,'color','w');
end
clear x
for ii=1
    
    f1=figure(15)
    
    
    dum=find(~isnan((arc.value(:,ii))));
    x=range(dum);
    
    y=arc.value(dum,ii);
    
    [fitresult] = fit( x, y, 'poly1' );
    h=plot(fitresult,x,y);
    
    hold on
    plot(x,y,'k.','MarkerSize',10,'HandleVisibility','off');
    legend('off')
    [r(1,ii),p(1,ii)]=corrcoef(x,y);
    
    %     legend(h,['r=',num2str(r)],'box','off')
    box('off')
    
    xlabel('IQ range')
    ylabel(sprintf( arc.label{ii,1}))
    
    
    set(gca,'FontSize',12)
    f1.Units ='centimeters';
    %     f1.OuterPosition= [10, 10, 8, 10];
    set(f1,'color','w');
    
    %     y2=plot(x(sig),y(sig),'ko','MarkerSize',8);
    
end
