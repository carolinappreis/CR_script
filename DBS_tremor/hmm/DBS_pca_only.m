clear all
close all
iii=[1];

for numb=1;
    %     :length(iii);
    
    %     load(strcat('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\0',num2str(iii(numb)),'_RS_PS.mat'))
    load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/0',num2str(iii(numb)),'_RS_PS.mat'));
    
    
    for  in2=1:3; % analysing the "main tremor axis"
        clearvars -except iii numb ttall ampall ph_stim LS tt1 in2 SmrData pca_ax1 pca_ax2 pca_ax3
        
        cd('/Users/Carolina/Documents/GitHub/CR_script/DBS_tremor')
        %         cd('C:\Users\creis\Documents\GitHub\CR_script\DBS_tremor')
        DBS_find_cond;
        
        
        for hh=1:2;
            clear handup Pxx F frange Pxxrange Fpeak tremor_or dummy envelope phase frequency tt
            
            
            handup=[];
            for i=1:length(start{hh,1})
                handup=[handup start{hh,1}(i):ending{hh,1}(i)]; %#ok<*AGROW>
            end
            handup=sort(handup,'ascend');
            
            [Pxx,F]=pwelch(tremor2(handup),samplerate,[],samplerate,samplerate);
            
            frange=F(3:10);
            Pxxrange=Pxx(3:10);
            
            Fpeak=frange(find(Pxxrange==max(Pxxrange)));
            
            if (Fpeak-2)>=1
                [b,a]=butter(2,[(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
            else
                [b,a]=butter(2,[(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
            end
            tremor_or=filtfilt(b,a,tremor2)*10*9.81/0.5;
            % tremor_or=zscore(tremor_or);
            dummy=hilbert(tremor_or);
            envelope=sqrt((real(dummy).^2)+(imag(dummy).^2));
            phase=angle(dummy);
            frequency=(smooth((1000/(2*pi))*diff(unwrap(angle(dummy))),500))';
            
            tremor=(data(3,:));% %score(:,1)';%
            ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
            ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
            
            
            [bb,aa]=butter(2,[0.8/(0.5*samplerate) ],'low'); %15
            % tremor_or=zscore(tremor_or);
            dummy=hilbert(tremor_or);
            envelope=sqrt((real(dummy).^2)+(imag(dummy).^2));
            phase=angle(dummy);
            frequency=(smooth((1000/(2*pi))*diff(unwrap(angle(dummy))),500))';
            
            
            tremor=(data(3,:));% %score(:,1)';%
            ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
            ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
            tremorx(1:size(ts1.data,3))=ts1.data;
            filt_x=filtfilt(bb,aa,tremorx);
            
            tremor=(data(6,:));% %score(:,1)';%
            ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
            ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
            tremory(1:size(ts1.data,3))=ts1.data;
            filt_y=filtfilt(bb,aa,tremory);
            
            tremor=(data(7,:));% %score(:,1)';%
            ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
            ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
            tremorz(1:size(ts1.data,3))=ts1.data;
            filt_z=filtfilt(bb,aa,tremorz);
            
            timeor=0:1/samplerate:(size(tremorx,2)-1)/samplerate;
            tremorxf=filtfilt(b,a,tremorx);
            tremoryf=filtfilt(b,a,tremory);
            tremorzf=filtfilt(b,a,tremorz);
            
            
            for j=1:length(start{hh,1})
                if (~isnan(start{hh,1}(j)))
                    x=[tremorxf(start{hh,1}(j):ending{hh,1}(j));tremoryf(start{hh,1}(j):ending{hh,1}(j));tremorzf(start{hh,1}(j):ending{hh,1}(j))];
                    [pc,score,latent,tsquare] = pca(x');
                    xxx(j,1:3)=pc(1:3,1);
                    ma(j)=(find(abs(xxx(j,1:3))==max(abs(xxx(j,1:3)))));
                end
            end
            
            tremor_or2=NaN(length(start{hh,1}),1);
            tremor_or3=NaN(length(start{hh,1}),1);
            
            for i=1:length(start{hh,1})
                if (~isnan(start{hh,1}(i))&& ma(i)==in2)
                    tremor_or3(i,1)=mean(envelope(start{hh,1}(i)-1000:start{hh,1}(i)));
                    tremor_or2(i,1)=(mean(envelope(ending{hh,1}(i)-1000:ending{hh,1}(i)))-mean(envelope(start{hh,1}(i)-1000:start{hh,1}(i))))/mean(envelope(start{hh,1}(i)-1000:start{hh,1}(i)));
                    xx{hh,1}(i)= xx{hh,1}(i);
                else
                    tremor_or2(i,1)=NaN;
                    tremor_or3(i,1)=NaN;
                    xx{hh,1}(i)= NaN;
                end
            end
            
            tt=NaN(20,12);
            yy=xx{hh,1}(:);
            
            for i=1:12
                tt(1:sum(yy==i),i)=tremor_or2(find(yy==i));
                tt(tt==0)=NaN;
            end
            arc(hh,numb,:)=nanmedian(tt);
        end
        if in2==1
            pca_ax1=arc;
            
        elseif in2==2
            pca_ax2=arc;
            
        elseif in2==3
            pca_ax3=arc;
            
        end
    end
    %cd('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data')
    %     cd('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data')
    %     clearvars -except pca_ax1 pca_ax2 pca_ax3
    %     %     save('pca')
end



clear all

load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/responsecurves_pcaonly.mat')
load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','aegean','stone','blushred');

cl(1,:)=blushred;
cl(2,:)=aegean;


for pp=1:2
    
    S1=repmat(squeeze(pca_ax1(pp,:,:))',1,3);
    S2=repmat(squeeze(pca_ax2(pp,:,:))',1,3);
    S3=repmat(squeeze(pca_ax3(pp,:,:))',1,3);
    
    for ii=1:size(S1,1)
        for i=size(pca_ax1,3)+1:size(pca_ax1,3)*2
            pca_1(pp,ii,i-12)=sum(S1(ii,(i-1:i+1)))./length(S1(ii,(i-1:i+1)));
            pca_2(pp,ii,i-12)=sum(S2(ii,(i-1:i+1)))./length(S2(ii,(i-1:i+1)));
            pca_3(pp,ii,i-12)= sum(S3(ii,(i-1:i+1)))./length(S3(ii,(i-1:i+1)));
        end
    end
end


for ms=1:2
    fig=[squeeze(pca_1(ms,1,:))';squeeze(pca_2(ms,1,:))';squeeze(pca_3(ms,1,:))'];
    
    %%------- POSTER IMAGES
    f=figure()
    
    for i=1:size(fig,1)
        subplot(1,size(fig,1),i)
        bar(0:30:330,fig(i,:),'LineWidth',0.5,'FaceColor',cl(ms,:),'EdgeColor',cl(ms,:))
        hold on
        yline(0,'LineWidth',1)
         ylim([-0.8 0.8])
%         yticks([ -0.5:0.25:0.5])
        box('off')
        ylabel('Change in tremor severity')
        xlabel('Stimulation phase (degrees)')
        title(['Motor',num2str(ms),'Axis',num2str(i)],'FontSize',12)
        
        
        f.Units = 'centimeters';
        % f.OuterPosition= [10, 10, 10, 10];
        set(gca,'XTickLabelRotation',45)
        set(gca,'FontSize',12)
        set(f,'color','w');
    end
    clear fig
end





for ms=1:2
    fig=[squeeze(pca_ax1(ms,1,:))';squeeze(pca_ax2(ms,1,:))';squeeze(pca_ax3(ms,1,:))'];
    
    %%------- POSTER IMAGES
    f=figure()
    
    for i=1:size(fig,1)
        subplot(1,size(fig,1),i)
        bar(fig(i,:),'LineWidth',0.5,'FaceColor',cl(ms,:),'EdgeColor',cl(ms,:))
        yline(0,'LineWidth',1)
         ylim([-0.8 0.8])
%         yticks([ -0.5:0.25:0.5])
        box('off')
        ylabel('Change in tremor severity')
        xlabel('Stimulation phase (degrees)')
        title(['Motor',num2str(ms),'Axis',num2str(i)],'FontSize',12)
        
        
        f.Units = 'centimeters';
        % f.OuterPosition= [10, 10, 10, 10];
        set(gca,'XTickLabelRotation',45)
        set(gca,'FontSize',12)
        set(f,'color','w');
    end
    clear fig
end