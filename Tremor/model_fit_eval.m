clear all
close all

% load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\cleaned_rc12_noaddon.mat')
% load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\newnonstim10.mat')
% load('C:\Users\creis\Documents\GitHub\CR_script\colour_pal.mat','blushred','squash');
% cl=blushred;
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/cleaned_rc12_noaddon.mat')
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/newnonstim10.mat')


cl=[0.5 0.5 0.5];

main=[1 1 3 1 3 3 3 3 1 1];
for pp=1:size(tt1,1)
    for kk=1:size(tt1,2)
        
        curve=nanmedian(tt1{pp,kk},1);
        SO=repmat(curve,1,3);
        %         SOUT=repmat((squeeze(nostimout(pp,kk,:)))',1,3);
        for i=size(tt1{pp,kk},2)+1:size(tt1{pp,kk},2)*2
            smooth_c(1,i-12)=sum(SO(1,(i-1:i+1)))./length(SO(1,(i-1:i+1)));
            %             smooth_ns(1,i-12)=sum(SOUT(1,(i-1:i+1)))./length(SOUT(1,(i-1:i+1)));
        end
        smoo_all(pp,kk,:)=smooth_c;  clear smooth_c
        %         nsmoo_all(pp,kk,:)=smooth_ns; clear smooth_ns
    end
    smoo_main(pp,:)=squeeze(smoo_all(pp,1,:));
    %     nsmoo_main(pp,:)=squeeze(nsmoo_all(pp,1,:));
    s_main(pp,:)=nanmedian(tt1{pp,1},1);
    %     ns_main(pp,:)=squeeze(nostimout(pp,1,:));
end


cl=[0.5 0.5 0.5];
close all
for i= 1:10
    subplot(2,5,i)
    x = 1:12;
    y = s_main(i,:);
    bar(x,y,'FaceColor',[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.5)
    box('off')
    hold on
    y=y';x=x';
    xData=x(find(~isnan(y))); yData=y(~isnan(y));
    mdp1=fitlm(x,y,'poly1')
    mdp2=fitlm(x,y,'poly2')
    mdp3=fitlm(x,y,'poly3')
    mdp4=fitlm(x,y,'poly4')
    mdk= fitlm(x,y,'constant');
    plot(mdp1.Fitted,'k','LineWidth',1.5)
    plot(mdp2.Fitted,'b','LineWidth',1.5)
    plot(mdp3.Fitted,'r','LineWidth',1.5)
    plot(mdp4.Fitted,'g','LineWidth',1.5)
    %plot(mdk.Fitted,'k','LineWidth',1.5)
    xticks([1:12])
    xticklabels([0:30:330])
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontSize',5)
    b = get(gca,'YTickLabel');
    set(gca,'YTickLabel',b,'FontSize',10)
    
    box('off')
    ylabel({'Change in tremor severity'},'FontSize',12)
    xlabel({'Stimulation phase (degrees)'},'FontSize',12)
    set(gca,'FontName','Arial','XTickLabelRotation',50)
    
    
    rsqr(i,:)=[ mdp1.Rsquared.Adjusted mdp2.Rsquared.Adjusted mdp3.Rsquared.Adjusted  mdp4.Rsquared.Adjusted ];
    F_stat(i,:)=[ coefTest(mdp1) coefTest(mdp2)  coefTest(mdp3)  coefTest(mdp4)];
    [r_val r_idx]=sort(rsqr(i,:),'descend');
    win_m(i,:)=[r_idx(1) r_val(1) F_stat(i,r_idx(1))];
    
    
    mcomp (i,:)=[mdk.ModelCriterion.AICc mdp1.ModelCriterion.AICc mdp2.ModelCriterion.AICc mdp3.ModelCriterion.AICc mdp4.ModelCriterion.AICc];
    [s_val s_idx]=sort(mcomp(i,:),'ascend');
    win(i,1)=s_idx(1);
    dif_aic(i,:)=[(abs(s_val(1)-s_val(2)))>2 (abs(s_val(1)-s_val(3)))>2 (abs(s_val(1)-s_val(4)))>2 (abs(s_val(1)-s_val(5)))>2];
end






%  beta0 = [(max(y)-min(y)) length(x)/2 1];
%     mg = @(F,x)(F(1)*exp(-((x-F(2))/F(3)).^2));
%     mdg = fitnlm(x,y,mg,beta0)

