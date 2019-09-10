%%% check start and end points and cf. with notes.

clear all
close all
iii=[1];

for numb=1;
    %     :length(iii);
    clearvars -except iii numb ttall ampall ph_stim LS tt1
          load(strcat('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\0',num2str(iii(numb)),'_RS_PS.mat'))
%     load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/0',num2str(iii(numb)),'_RS_PS.mat'));

    in2=1; % analysing the "main tremor axis"
    
    DBS_find_cond;
    clear handup Pxx F frange Pxxrange Fpeak tremor_or dummy envelope phase frequency
    
    
    for hh=1:2;
        
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
        
        
        [b,a]=butter(2,[0.8/(0.5*samplerate) ],'low'); %15
        % tremor_or=zscore(tremor_or);
        dummy=hilbert(tremor_or);
        envelope=sqrt((real(dummy).^2)+(imag(dummy).^2));
        phase=angle(dummy);
        frequency=(smooth((1000/(2*pi))*diff(unwrap(angle(dummy))),500))';
        tremor=(data(3,:));% %score(:,1)';%
        ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
        ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
        tremorx(1:size(ts1.data,3))=ts1.data;
        filt_x=filtfilt(b,a,tremorx);
        tremor=(data(6,:));% %score(:,1)';%
        ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
        ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
        tremory(1:size(ts1.data,3))=ts1.data;
        filt_y=filtfilt(b,a,tremory);
        tremor=(data(7,:));% %score(:,1)';%
        ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
        ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
        tremorz(1:size(ts1.data,3))=ts1.data;
        filt_z=filtfilt(b,a,tremorz);
        timeor=0:1/samplerate:(size(tremorx,2)-1)/samplerate;
        
        %         subplot(3,1,1)
        %         plot(timeor(1,:), filt_x(1,:));
        %         subplot(3,1,2)
        %         plot(timeor(1,:), filt_y(1,:));
        %         subplot(3,1,3)
        %         plot(timeor(1,:), filt_z(1,:));
        %
        %     subplot(3,1,1)
        %     plot(timeor(1,:), tremorx(1,:));
        %     subplot(3,1,2)
        %     plot(timeor(1,:), tremory(1,:));
        %     subplot(3,1,3)
        %     plot(timeor(1,:), tremorz(1,:));
        
%         dt1_s=floor((dt1_s./samplerateold)*samplerate)+addon;
%         dt1_e=floor((dt1_e./samplerateold)*samplerate)+addon+addon_end;%floor(5*samplerate);
%         dt2_s=floor((dt2_s./samplerateold)*samplerate)+addon;
%         dt2_e=floor((dt2_e./samplerateold)*samplerate)+addon+addon_end;%floor(5*samplerate);
%         
%--------------check         
%         for i=1:length(dt1_e)
%             figure(1)
%             subplot(2,1,1)
%             plot3(timeor(1,dt1_s(i):dt1_e(i)),filt_y(1,dt1_s(i):dt1_e(i)), filt_z(1,dt1_s(i):dt1_e(i)));
%             hold on
%             subplot(2,1,2)
%             plot3(timeor(1,dt2_s(i):dt2_e(i)),filt_y(1,dt2_s(i):dt2_e(i)), filt_z(1,dt2_s(i):dt2_e(i)));
%             hold on
%         end
%         
        
        
        tremor_or2=NaN(length(start{hh,1}),1);
        tremor_or3=NaN(length(start{hh,1}),1);
        
        for i=1:length(start{hh,1})
            if (~isnan(start{hh,1}(i)))
                tremor_or3(i,1)=mean(envelope(start{hh,1}(i)-1000:start{hh,1}(i)));
                tremor_or2(i,1)=(mean(envelope(ending{hh,1}(i)-1000:ending{hh,1}(i)))-mean(envelope(start{hh,1}(i)-1000:start{hh,1}(i))))/mean(envelope(start{hh,1}(i)-1000:start{hh,1}(i)));
                xx{hh,1}(i)= xx{hh,1}(i);
            else
                tremor_or2(i,1)=NaN;
                tremor_or3(i,1)=NaN;
                xx{hh,1}(i)= NaN;
            end
        end
        
        %         %% criteria for outliers
        %
        %         idx_outl=find(tremor_or2>(nanmean(tremor_or2)+(2*(nanstd(tremor_or2))))|tremor_or2<(nanmean(tremor_or2)-(2*(nanstd(tremor_or2)))));
        %         tremor_or2(idx_outl,1)=NaN;
        %         tremor_or3(idx_outl,1)=NaN;
        %         xx(1,idx_outl)=NaN;
        
        
        tt=NaN(25,12);
        amp=NaN(25,12);
        yy=xx{hh,1}(:);
        
        for i=1:12
            tt(1:sum(yy==i),i)=tremor_or2(find(yy==i));
            amp(1:sum(yy==i),i)=tremor_or3(find(yy==i));
            ph_stim{hh,1}(numb,i)=sum(yy==i);
        end
        
        clear yy;
        
        tt1{hh,1}{numb,1}=tt;
        
        ttall {hh,1}(numb,:)=nanmedian(tt);
        ampall {hh,1}(numb,:)=nanmedian(amp);
        
        % for s=1:size(tt,2)
        %     for i =1:100000;
        %         yy1=xx(randperm(size(xx,2)));
        %         tt2(1:sum(yy1==s),1)=tremor_or2(find(yy1==s));
        %         tt3(i,s)=nanmedian(tt2,1);
        %         clear tt2
        %     end
        % end
        % LS (numb,:,:)=tt3;
        
        for rr=1:100000
            LS{hh,1}(numb,rr)=nanmedian(tt(randi(length(start{hh,1}),1,10)));
        end
 
        
%--------------check        
%         figure(hh+2)
%         for f=1:6
%             epoch=start{hh,1}(f):ending{hh,1}(f);
%             plot3(filt_x(1,epoch),filt_y(1,epoch), filt_z(1,epoch),'color',rand(1,3));
%             hold on
%         end
%         xlabel('x axis')
%         ylabel('y axis')
%         zlabel('z axis')
        
    end
    clearvars -except ttall iii numb ampall ph_stim LS tt1 hh
    
end
clearvars -except ttall ampall ph_stim LS tt1
cd('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data')
% cd('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data')
% save('DBS_amp_ARC.mat')

% cd('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data')
load('DBS_amp_ARC.mat')

load('C:\Users\creis\Documents\GitHub\CR_script\colour_pal.mat','blushred','squash');
% load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','squash')
cl=blushred;
cl1=squash;

figure()
subplot(2,1,1)
bar(ttall{1,1})
hold on
plot((cell2mat(tt1{1,1}(:)))','.')
box('off')
title('posture')
subplot(2,1,2)
bar(ttall{2,1})
hold on
plot((cell2mat(tt1{2,1}(:)))','.')
box('off')
title('spiral')


%%%  smooth
P=repmat(ttall{1,1},1,3);
S=repmat(ttall{2,1},1,3);


for ii=1:size(ttall{1,1},1)
    for i=size(ttall{1,1},2)+1:size(ttall{1,1},2)*2
        pst_s(ii,i-12)=sum(P(ii,(i-1:i+1)))./length(P(ii,(i-1:i+1)));
        sprl_s(ii,i-12)=sum(S(ii,(i-1:i+1)))./length(S(ii,(i-1:i+1)));
    end
end

figure()
subplot(2,1,1)
bar(pst_s,'FaceColor',cl,'EdgeColor',cl)
box('off')
title('posture')
subplot(2,1,2)
bar(sprl_s,'FaceColor',cl,'EdgeColor',cl)
box('off')
title('spiral')


y=[pst_s;sprl_s];
for i =1:2;
    clearvars -except i y rs_gauss rs_sin cl smo_s
    figure(i)
    subplot(1,2,1)
    bar(y(i,:),'FaceColor',cl,'EdgeColor',cl)
    hold on
    rsg=gauss_fit(y(i,:));
    rs_gauss(i,:)=rsg.adjrsquare;
    %     legend( 'ARC', 'gaussian fit', 'Location', 'NorthEast', 'Interpreter', 'none');
    %     legend('boxoff')
    %    xlabel( 'Stim phase', 'Interpreter', 'none' );
    %    ylabel( 'Amplitude change', 'Interpreter', 'none' );
    title(['adjr^2: ',num2str(rs_gauss(i))])
    %   ylim([-(max(abs(y)))-0.05 (max(abs(y))+0.05)])
    xlabel('');ylabel('');legend('off')
    box('off')
    
    subplot(1,2,2)
    bar(y(i,:),'FaceColor',cl,'EdgeColor',cl)
    hold on
    rss=sin_fit(y(i,:));
    rs_sin(i,:)=rss.adjrsquare;
    xlabel('');ylabel('');legend('off')
    %     legend( 'ARC','Sin fit', 'Location', 'NorthEast', 'Interpreter', 'none');
    %     legend('boxoff')
    %     xlabel( 'Stim phase', 'Interpreter', 'none' );
    %     ylabel( 'Amplitude change', 'Interpreter', 'none' );
    title(['adjr^2: ',num2str(rs_sin(i))])
    %     ylim([-(max(abs(y)))-0.05 (max(abs(y))+0.05)])
    box('off')
    hold off
end

