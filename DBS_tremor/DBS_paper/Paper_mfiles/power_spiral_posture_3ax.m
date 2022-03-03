cd('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA')
cd('/Users/Carolina/Documents/GitHub/CR_script/DBS_tremor/DBS_paper')
clear; close
cohort = [ 1 3 4 6];
cond={'NS';'RS'};
clust=struct; out=struct; start=cell(10,1); ending=cell(10,1); yy=cell(10,3); h_up=cell(10,3); s=struct; freq_bl=[]; amp_bl=[];
rng('default')
gen=(rng);
names={'Z';'Y';'X'};
type=0;


for iii =  1:length(cohort)
    clearvars -except ml type names cohort cond iii clust s start ending yy out gen h_up spiral freq_bl amp_bl avg power dat
    
    co=1;
    load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA/P0',num2str(cohort(iii)),'_',num2str(cond{co,1}),'.mat'))
    leg=[];
    for mmmm=1:2
        
        spiral=mmmm-1;
        
        [d]=dbs_preprocess(SmrData); samplerateold=d.samplerateold; samplerate=d.samplerate;
        [peak_ax, start, ending, yy, h_up]= srtend(d,samplerateold,samplerate,iii,co,start,ending,yy,h_up,spiral);
        [s]=zfiltenv_simple(d,peak_ax,co,iii,s,samplerate);
        [match_ax]=link_ax(spiral);
        
        f1=figure(1)
        subplot(2,2,iii)
        load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','aegean','stone','squash','sapphire','azure','ax2','ax3');
        color_b1=[stone;ax2;ax3];
        for aa=1:3
        [Pxx,F] = pwelch(s.raw{iii, 1} (aa,h_up{iii,1}), samplerate, [], samplerate, samplerate);
        if mmmm==1
        plot(F(3:10),Pxx(3:10),'-','Color',color_b1(aa,:),'LineWidth',1)
        else
        plot(F(3:10),Pxx(3:10),'--','Color',color_b1(aa,:),'LineWidth',1)
        end
        hold on
        box('off')
        title(sprintf('Patient %d',(iii)))
        end

%         signal=sqrt(((s.filth{iii, 1}(1,:)).^2)+((s.filth{iii, 1}(2,:)).^2)+((s.filth{iii, 1}(3,:)).^2));
%         %
%         if type==0
%             [Pxx,F] = pwelch(s.raw{iii, 1} (match_ax(co,iii,1),h_up{iii,1}), samplerate, [], samplerate, samplerate);
%         else
%             [Pxx,F] = pwelch(signal(h_up{iii,1}), samplerate, [], floor(samplerate*2), samplerate);
%         end
%         
%         power(iii,:)=Pxx;
%         f1=figure(1)
%         subplot(2,2,iii)
%         cut1=find(F==2);
%         cut2=find(F==10);
%         plot((F(cut1:cut2)),(power(iii,cut1:cut2)),'LineWidth',2) 
%         hold on
%         box('off')
%         title(sprintf('patient %d',(iii)))
%         legend('Posture','Spiral')
%         leg=[leg match_ax(co,iii,1)];
        
    end
%     ml(iii,:)=names(leg);

     legend ('P_Z','P_Y','P_X','S_Z','S_Y','S_X')
     legend('boxoff')
%     %         plot(log10(F),log10(power(iii,:)))
%     %         xlim([log10(2) 1])
%     %         xticks([log10(2):0.2:1])
%     %         xticklabels([  10^(0.2) 10^(0.4) 10^(0.6) 10^(0.8) 10])
    hold on
    set(gca,'FontSize',14)
    f1.Units ='centimeters';
    f1.OuterPosition= [10, 10, 20, 20];
    set(f1,'color','w');
    
end
