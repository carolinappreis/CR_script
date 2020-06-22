%%% check start and end points and cf. with notes.

clear all
close all
iii=[1];

for numb=1;
    %     :length(iii);
    clearvars -except iii numb ttall amp_1b amp_1l ph_stim LS tt1
%            load(strcat('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\0',num2str(iii(numb)),'_RS_PS.mat'))

   for  in2=1:3; % analysing the "main tremor axis"
    load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/0',num2str(iii(numb)),'_RS_PS.mat'));

    
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
        
        
        tremor_or2=NaN(length(start{hh,1}),1);
        amp_bef=NaN(length(start{hh,1}),1);
        amp_last=NaN(length(start{hh,1}),1);
        
        
        for i=1:length(start{hh,1})
            if (~isnan(start{hh,1}(i)))
                amp_bef(i,1)=mean(envelope(start{hh,1}(i)-1000:start{hh,1}(i)));
                amp_last(i,1)=mean(envelope(ending{hh,1}(i)-1000:ending{hh,1}(i)));
                tremor_or2(i,1)=(mean(envelope(ending{hh,1}(i)-1000:ending{hh,1}(i)))-mean(envelope(start{hh,1}(i)-1000:start{hh,1}(i))))/mean(envelope(start{hh,1}(i)-1000:start{hh,1}(i)));
%                 plot(envelope(start{hh,1}(i)-1000:ending{hh,1}(i)));
                xx{hh,1}(i)= xx{hh,1}(i);
            else
                tremor_or2(i,1)=NaN;
                amp_bef(i,1)=NaN;
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
        amp_b=NaN(25,12);
        amp_l=NaN(25,12);
        yy=xx{hh,1}(:);
        
        for i=1:12
            tt(1:sum(yy==i),i)=tremor_or2(find(yy==i));
            amp_b(1:sum(yy==i),i)=amp_bef(find(yy==i));
            amp_l(1:sum(yy==i),i)=amp_last(find(yy==i));
            ph_stim{hh,1}(in2,i)=sum(yy==i);
        end
        
        clear yy;
        
        tt1{in2,1}{hh,1}=tt;
        
        ttall (in2,hh,:)=nanmedian(tt);
        amp_1b (in2,hh,:)=nanmedian(amp_b);
        amp_1l(in2,hh,:)=nanmedian(amp_l);       
        
        for rr=1:100000
            LS(in2,hh,rr)=nanmedian(tt(randi(length(start{hh,1}),1,10)));
        end
        
    end
           clearvars -except ttall iii numb amp_1b amp_1l ph_stim LS tt1 hh in2

    
   end

   
end

clear numb in2 iii hh 
cd('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data')
% save('3ax_all')
clearvars -except ttall ampall ph_stim LS tt1

%  load('C:\Users\creis\Documents\GitHub\CR_script\colour_pal.mat','blushred','aegean');
 load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','aegean')
cl=aegean;
cl1=blushred;

figure()
subplot(2,1,1)
bar(ttall{1,1},'FaceColor',cl,'EdgeColor',cl)
hold on
% plot((cell2mat(tt1{1,1}(:)))','.')
xlabel('stimulation phase')
ylabel ('change tremor severity')
box('off')
title('posture')
subplot(2,1,2)
bar(ttall{2,1},'FaceColor',cl1,'EdgeColor',cl1)
hold on
% plot((cell2mat(tt1{2,1}(:)))','.')
xlabel('stimulation phase')
ylabel ('change tremor severity')
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
bar(0:30:330,pst_s,'FaceColor',cl,'EdgeColor',cl)
ax = gca; ax.FontSize = 12; ax.YLim = [-0.4 0.4];
box('off')
xlabel('stimulation phase','FontSize',14)
ylabel ('change tremor severity','FontSize',14)
% title('Posture','FontSize',14)
subplot(2,1,2)
bar(0:30:330,sprl_s,'FaceColor',cl1,'EdgeColor',cl1)
ax = gca; ax.FontSize = 12; ax.YLim = [-0.4 0.4];
xlabel('stimulation phase','FontSize',14)
ylabel ('change tremor severity','FontSize',14)
box('off')
% title('Spiral', 'FontSize',14)



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

