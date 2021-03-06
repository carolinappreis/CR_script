%%% check start and end points and cf. with notes.

clear all
% close all
iii=[1];

for numb=1;
    %     :length(iii);
    clearvars -except iii numb ttall ampall ph_stim LS tt1
%      load(strcat('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\0',num2str(iii(numb)),'_RS_PS.mat'))
    load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/0',num2str(iii(numb)),'_RS_PS.mat'));


    in2=1; % analysing the "main tremor axis"

    if in2==1
        in=3;
    elseif in2==2 % other axis 1
        in=6;
    elseif in2==3 % other axis 2
        in=7;
    end
    data=SmrData.WvData;
    samplerateold=SmrData.SR;
    tremor=(data(in,:));
    addon=92; addon_end=35;

    time=0:1/samplerateold:(size(data,2)-1)/samplerateold;

    %%% downsample

    ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
    ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
    tremor2(1:size(ts1.data,3))=ts1.data;
    samplerate=1000;

    %%% determine stimulation time points
    index=[];
    for i=2:size(data,2)-1
        if data(2,i-1)<2.5 && data(2,i)>2.5
            index=[index i];
        end
    end
    clear i

    indexes4=[index(1) index(find(diff(index)./samplerateold > 0.95)+1)];
    indexes3=[index(find(diff(index)./samplerateold > 0.95)) index(end)];

    dd2=round(data(4,:)*100)./100;
    for i=1:length(indexes4)
        xx(i)=round(dd2(indexes4(i))./0.1); %#ok<*SAGROW>
    end
    clear i

    start=floor((indexes4./samplerateold)*samplerate)+addon;
    ending=floor((indexes3./samplerateold)*samplerate)+addon+addon_end;%floor(5*samplerate);

    %%% when patient's hand is up
    handup=[];
    for i=1:length(start)
        handup=[handup start(i):ending(i)]; %#ok<*AGROW>
    end
    clear i
    handup=sort(handup,'ascend');


    %%% tremor characteristics
    [Pxx,F]=pwelch(tremor2(handup),samplerate,[],samplerate,samplerate);

    frange=F(3:10);
    Pxxrange=Pxx(3:10);

    Fpeak=frange(find(Pxxrange==max(Pxxrange))); %#ok<*FNDSB>

    if (Fpeak-2)>=1
        [b,a]=butter(2,[(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
    else
        [b,a]=butter(2,[(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
    end

    tremor_or=filtfilt(b,a,tremor2)*10*9.81/0.5;

    dummy=hilbert(tremor_or);
    envelope=sqrt((real(dummy).^2)+(imag(dummy).^2));

    th1=(Fpeak*5)./2;
    th2=(Fpeak*5)+round((Fpeak*5)./5);

    %-------

    new=find(data(2,:)>4);
    difp=find((diff(new))>100000);
    ep_1=[new(difp) new(end)];
    sp_1=[new(1) new(difp+1)];


    if numb==1
        dt1_s=[sp_1(2:2:end)];dt1_e=[ep_1(2:2:end)];
        dt2_s=sp_1(3:2:end); dt2_e=ep_1(3:2:end);
        sp=sp_1(1,2:end);
        ep=ep_1(1,2:end);
    end
    %
    %         plot(time,data(4,:))
    %     hold on
    %     plot(time(sp),data(4,sp),'r.')
    %     plot(time(ep),data(4,ep),'k.')


    for ik=1:length(sp) %%find double start and end points in a stimulation run

        s=(find(([indexes4>=sp(ik)]+[indexes4<=ep(ik)])==2));
        e=(find(([indexes3>=sp(ik)]+[indexes3<=ep(ik)])==2));
        tks=(find(diff(xx(s))==0))+1;
        tke=(find(diff(xx(e))==0));

        indexes4(s(tks))=NaN;
        indexes3(e(tke))=NaN;
        xx(e(tke))=NaN;

    end

    %%%%% find runs with trigering issues (too few, too many pulses)
    %         for it=1:length(indexes4)
    %             if numel(index(find(index==indexes4(it)):find(index==indexes3(it))))>=th1 && numel(index(find(index==indexes4(it)):find(index==indexes3(it))))<=th2
    %                 indexes4(it)=indexes4(it);
    %                 indexes3(it)=indexes3(it);
    %                 xx(it)=xx(it);
    %             else
    %                 indexes4(it)=NaN;
    %                 indexes3(it)=NaN;
    %                 xx(it)=NaN;
    %             end
    %         end
    %%%%%%%%%%%%%%%%
    indexes4=indexes4(~isnan(indexes4));
    indexes3=indexes3(~isnan(indexes3));
    xx=xx(~isnan(xx));


    start1=[];
    ending1=[];
    xx1=[];
    for il=1:length(dt1_s)
        start1=[start1 indexes4(find(([indexes4>=dt1_s(il)]+[indexes4<=dt1_e(il)])==2))];
        ending1=[ending1 indexes3(find(([indexes3>=dt1_s(il)]+[indexes3<=dt1_e(il)])==2))];
        xx1=[xx1 xx(find(([indexes4>=dt1_s(il)]+[indexes4<=dt1_e(il)])==2))];
    end

    start2=[];
    ending2=[];
    xx2=[];
    for il=1:length(dt2_s)
        if numb==1&&il==3 | numb==1&&il==9
            dums=indexes4(find(([indexes4>=dt2_s(il)]+[indexes4<=dt2_e(il)])==2));
            start2=[start2 dums(1,3:end)]; clear dums
            dume=indexes3(find(([indexes3>=dt2_s(il)]+[indexes3<=dt2_e(il)])==2));
            ending2=[ending2 dume(1,3:end)]; clear dume
            dumx=xx(find(([indexes4>=dt2_s(il)]+[indexes4<=dt2_e(il)])==2));
            xx2=[xx2 dumx(1,3:end)];clear dumx
        else
            dums=indexes4(find(([indexes4>=dt2_s(il)]+[indexes4<=dt2_e(il)])==2));
            start2=[start2 dums]; clear dums
            dume=indexes3(find(([indexes3>=dt2_s(il)]+[indexes3<=dt2_e(il)])==2));
            ending2=[ending2 dume]; clear dume
            dumx=xx(find(([indexes4>=dt2_s(il)]+[indexes4<=dt2_e(il)])==2));
            xx2=[xx2 dumx];clear dumx

        end
    end
    %
    %     figure()
    %     plot(time,data(4,:))
    %     hold on
    %     plot(time(index),data(4,index),'r.')
    %     plot(time(start1),data(4,start1),'ko')
    %     plot(time(ending1),data(4,ending1),'bo')


    clear start ending
    start{1,1}=floor((start1./samplerateold)*samplerate)+addon;
    ending{1,1}=floor((ending1./samplerateold)*samplerate)+addon+addon_end;%floor(5*samplerate);
    start{2,1}=floor((start2./samplerateold)*samplerate)+addon;
    ending{2,1}=floor((ending2./samplerateold)*samplerate)+addon+addon_end;%floor(5*samplerate);
    clear xx
    xx{1,1}=xx1;
    xx{2,1}=xx2;

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

        dt1_s=floor((dt1_s./samplerateold)*samplerate)+addon;
        dt1_e=floor((dt1_e./samplerateold)*samplerate)+addon+addon_end;%floor(5*samplerate);
        dt2_s=floor((dt2_s./samplerateold)*samplerate)+addon;
        dt2_e=floor((dt2_e./samplerateold)*samplerate)+addon+addon_end;%floor(5*samplerate);


        %         for i=1:length(dt1_e)
        %             figure(1)
        %             subplot(2,1,1)
        %             plot3(timeor(1,dt1_s(i):dt1_e(i)),tremorx(1,dt1_s(i):dt1_e(i)), tremorz(1,dt1_s(i):dt1_e(i)));
        %             hold on
        %             subplot(2,1,2)
        %             plot3(timeor(1,dt2_s(i):dt2_e(i)),tremorx(1,dt2_s(i):dt2_e(i)), tremorz(1,dt2_s(i):dt2_e(i)));
        %             hold on
        %         end
        %


        tremor_or2=NaN(length(start{hh,1}),1);
        tremor_or3=NaN(length(start{hh,1}),1);

        % amplitude
        tremor_or2=NaN(length(start{hh,1}),1);
        for i=1:length(start{hh,1})
            if (~isnan(start{hh,1}(i)))
                tremor_or3(i,1)=mean(envelope(start{hh,1}(i)-1000:start{hh,1}(i)));
                tremor_or2(i,1)=(mean(envelope(ending{hh,1}(i)-1000:ending{hh,1}(i)))-mean(envelope(start{hh,1}(i)-1000:start{hh,1}(i))))/mean(envelope(start{hh,1}(i)-1000:start{hh,1}(i)));

            else
                tremor_or2(i,1)=NaN;
                tremor_or3(i,1)=NaN;
            end
        end


        amp_1=NaN(2,round(size(tremor_or3,1)./2));
        ch_a1=NaN(2,round(size(tremor_or3,1)./2));
        pha_idx=NaN(2,round(size(tremor_or3,1)./2));
        m=1;
        n=1;
        for i=1:length(start{hh,1})
            if tremor_or3(i)<=nanmedian(tremor_or3)
                amp_1(1,n)= tremor_or3(i);
                ch_a1(1,n)= tremor_or2(i);
                pha_idx(1,n)=xx{hh,1}(i);
                n=n+1;

            else
                amp_1(2,m)= tremor_or3(i);
                ch_a1(2,m)= tremor_or2(i);
                pha_idx(2,m)=xx{hh,1}(i);
                m=m+1;

            end
        end

        clear tt1 tt2 amp1 amp2
        tt1=NaN(20,12);
        amp1=NaN(20,12);
        tt2=NaN(20,12);
        amp2=NaN(20,12);

        for i=1:12
            tt1(1:sum(pha_idx(1,:)==i),i)=ch_a1(1,find(pha_idx(1,:)==i));
            amp1(1:sum(pha_idx(1,:)==i),i)=amp_1(1,find(pha_idx(1,:)==i));
            tt2(1:sum(pha_idx(2,:)==i),i)=ch_a1(2,find(pha_idx(2,:)==i));
            amp2(1:sum(pha_idx(2,:)==i),i)=amp_1(2,find(pha_idx(2,:)==i));
        end

        effect1=nanmedian(tt1);
        effect2=nanmedian(tt2);

        ef_1=repmat(effect1,1,3);
        ef_2=repmat(effect2,1,3);

        for i=size(effect1,2)+1:size(effect1,2)*2
            arc1{hh,1}(numb,i-12)=nansum(ef_1(1,(i-1:i+1)))./length(ef_1(1,(i-1:i+1)));
            arc2{hh,1}(numb,i-12)=nansum(ef_2(1,(i-1:i+1)))./length(ef_2(1,(i-1:i+1)));
        end
        %
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

end
 clearvars -except arc1 arc2
% % cd('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data')
%  cd('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA')
% 
%  save('arc_mediansplit.mat')
%

%%________________________________________________________________________
% clear all
% close all
% cd('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data')
% cd('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data')
% 
% load('arc_mediansplit.mat')

 load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','aegean','grey')
%  load('C:\Users\creis\Documents\GitHub\CR_script\colour_pal.mat','blushred','aegean','grey');


 cl2=[0.7 0.7 0.7];


% 
% for hh=1:size(arc1,1)
%     subplot(1,size(arc1,1),hh)
%     bar(0:30:330,arc2{hh,1},'FaceColor',cl1,'EdgeColor',cl1)
%     hold on
%     bar(0:30:330,arc1{hh,1},'LineStyle','--','LineWidth',1,'FaceColor','none','EdgeColor',cl2)
%     ax = gca; ax.FontSize = 12; ax.YLim = [-0.8 0.8];
%     xlabel('stimulation phase','FontSize',14)
%     ylabel ('change tremor severity','FontSize',14)
%     box('off')
%     if hh==1
%         title('posture')
%     else
%         title('spiral')
%     end
% end


hh=2;
%  cl1=blushred;
%   cl1=aegean;

% f=figure()
% bar(0:30:330,arc2{hh,1},'LineWidth',2,'FaceColor',cl1,'EdgeColor',cl1)
% hold on
% bar(0:30:330,arc1{hh,1},'LineWidth',2,'FaceColor','none','EdgeColor',cl2)
% yline(0,'LineWidth',1)
% ylim([-0.75 0.75])
% yticks([ -1:0.25:1])
% box('off')
% ylabel('Change in tremor severity')
% xlabel('Stimulation phase (degrees)')
%     
% f.Units = 'centimeters';
% f.OuterPosition= [10, 10, 12, 12];
% set(gca,'XTickLabelRotation',45)
% set(gca,'FontSize',14)
% set(f,'color','w');

f=figure()
y=[arc1{hh,1};arc2{hh,1}]';
b = bar(0:30:330,y);
yline(0,'LineWidth',1)
ylim([-0.75 0.75])
yticks([ -1:0.25:1])
box('off')
ylabel('Change in tremor severity')
xlabel('Stimulation phase (degrees)')
    
f.Units = 'centimeters';
f.OuterPosition= [10, 10, 12, 12];
set(gca,'XTickLabelRotation',45)
set(gca,'FontSize',14)
set(f,'color','w');



