function   [out]=s_preprocess(SmrData, out, iii, co)

data=SmrData.WvData;
out.samplerate=SmrData.SR;
data_t=data([3 6 7],:);
out.samplerateold=1000;

if co==1
       
    NS_BE_S
    out.b_trials{iii,co}=floor((hu{iii,:}./out.samplerateold)*out.samplerate);
    out.e_trials{iii,co}=floor((hd{iii,:}./out.samplerateold)*out.samplerate);
    
    yy{iii,co}=NaN;
    
    %----------------------------------------------------------------------------------------------
else
    addon=92; addon_end=35;
    
    % determine stimulation time points
    index = [];
    for i = 2:size(data, 2)-1
        if data(2,i-1)<2.5 && data(2,i)>2.5
            index = [index i];
        end
    end
    clear i
    
    % Find trigger
    indexes4 = [index(1) index(find(diff(index) ./ out.samplerateold > 0.95)+1)];
    indexes3 = [index(find(diff(index) ./ out.samplerateold > 0.95)) index(end)];
    
    
    dd2 = round(data(4, :) *100) ./ 100;
    for i = 1:length(indexes4)
        xx(i) = round(dd2(indexes4(i)) ./ 0.1); %#ok<*SAGROW>
    end
    clear i
    
    
    new = find(data(2,:) > 4);
    difp = find((diff(new)) > 100000); % are you trying to threshold at 9.6 seconds?
    ep_1 = [new(difp) new(end)];
    sp_1 = [new(1) new(difp+1)];
    
    %%% input start all trial
    start_t=1;
    
    sp=sp_1(1,start_t:end);
    ep=ep_1(1,start_t:end);
    
    %%% posture/spiral trials
    dt1_s=[sp_1(start_t:2:end)];dt1_e=[ep_1(start_t:2:end)];
    dt2_s=sp_1(start_t+1:2:end); dt2_e=ep_1(start_t+1:2:end);
    
    %     time=1:length(data(1,:));
    %     plot(time,data(4,:))
    %     hold on
    %     plot(time(dt2_s),data(4,dt2_s),'r.')
    %     plot(time(dt2_e),data(4,dt2_e),'k.')
    
    out.b_trials{iii,co}=dt2_s+floor((addon./out.samplerateold)*out.samplerate);
    out.e_trials{iii,co}=dt2_e+floor(((addon+addon_end)./out.samplerateold)*out.samplerate);
    
    
    for ik=1:length(sp) %%find double start and end points in a stimulation run
        
        s=(find(([indexes4>=sp(ik)]+[indexes4<=ep(ik)])==2));
        e=(find(([indexes3>=sp(ik)]+[indexes3<=ep(ik)])==2));
        tks=(find(diff(xx(s))==0))+1;
        tke=(find(diff(xx(e))==0));
        
        indexes4(s(tks))=NaN;
        indexes3(e(tke))=NaN;
        xx(e(tke))=NaN;
        
    end
    
    
    %%%% find runs with trigering issues (too few, too many pulses)
    % % %     th1=(Fpeak*5)./2;
    % % %     th2=(Fpeak*5)+round((Fpeak*5)./5);
    % % %     for it=1:length(indexes4)
    % % %         if numel(index(find(index==indexes4(it)):find(index==indexes3(it))))>=th1 && numel(index(find(index==indexes4(it)):find(index==indexes3(it))))<=th2
    % % %             indexes4(it)=indexes4(it);
    % % %             indexes3(it)=indexes3(it);
    % % %             xx(it)=xx(it);
    % % %         else
    % % %             indexes4(it)=NaN;
    % % %             indexes3(it)=NaN;
    % % %             xx(it)=NaN;
    % % %         end
    % % %     end
    % % %     %%%%%%%%%%%%%%%
    % % %     indexes4=indexes4(~isnan(indexes4));
    % % %     indexes3=indexes3(~isnan(indexes3));
    % % %     xx=xx(~isnan(xx));
    
    
    start1=[];
    ending1=[];
    xx1=[];
    for il=1:length(dt1_s)
        start1=[start1 indexes4(find(([indexes4>=dt1_s(il)]+[indexes4<=dt1_e(il)])==2))]; % intersect([1 2 3],[3 4 5])
        ending1=[ending1 indexes3(find(([indexes3>=dt1_s(il)]+[indexes3<=dt1_e(il)])==2))];
        xx1=[xx1 xx(find(([indexes4>=dt1_s(il)]+[indexes4<=dt1_e(il)])==2))];
    end
    
    start2=[];
    ending2=[];
    xx2=[];
    for il=1:length(dt2_s)
        dums=indexes4(find(([indexes4>=dt2_s(il)]+[indexes4<=dt2_e(il)])==2));
        start2=[start2 dums]; clear dums
        dume=indexes3(find(([indexes3>=dt2_s(il)]+[indexes3<=dt2_e(il)])==2));
        ending2=[ending2 dume]; clear dume
        dumx=xx(find(([indexes4>=dt2_s(il)]+[indexes4<=dt2_e(il)])==2));
        xx2=[xx2 dumx];clear dumx
    end
    
    
% %         figure()
% %         time=1:length(data(4,:));
% %         plot(time,data(4,:))
% %         hold on
% %         plot(time(index),data(4,index),'r.')
% %         plot(time(start2),data(4,start2),'ko')
% %         plot(time(ending2),data(4,ending2),'bo')
    
    
    clear start ending xx
    pstart{1,1}=start1+floor((addon./out.samplerateold)*out.samplerate);
    pending{1,1}=ending1+floor(((addon+addon_end)./out.samplerateold)*out.samplerate);
    pstart{2,1}=start2+floor((addon./out.samplerateold)*out.samplerate);
    pending{2,1}=ending2+floor(((addon+addon_end)./out.samplerateold)*out.samplerate);
    
    xx{1,1}=xx1;
    xx{2,1}=xx2;
    
    %%% choosing start{1,1}/ending{1,1}/xx{1,1} to get posture only _ check!
    
    out.start{iii,co}= pstart{2,1};
    out.ending{iii,co} = pending{2,1};
    yy{iii,co} = xx{2,1};
    
% %     figure()
% %     index_ds=index;
% %     time=1:length(data_t);
% %     plot(time,data_t(1,:))
% %     hold on
% %     plot(time(1,index_ds),data_t(1,index_ds),'r.')
% %     plot(time(start{iii,co}),data_t(1,start{iii,co}),'bo')
% %     plot(time(ending{iii,co}),data_t(1,ending{iii,co}),'ko')
%     ----------------------------------------------------------------------------------------------
end

test=[];
handup = [];
for i = 1:length(out.b_trials{iii,co})
    handup = [handup out.b_trials{iii,co}(i):out.e_trials{iii,co}(i)]; %#ok<*AGROW>
end
clear i
out.test=sort(test,'ascend');
handup = sort(handup,'ascend');
out.h_up{iii,co}=handup;

for aa = 1:3
    [Pxx,F] = pwelch(data_t(aa,handup), out.samplerate, [], round(out.samplerate), out.samplerate);
    frange = F(3:10);
    Pxxrange = Pxx(3:10);
    Freqpeak(aa,:) = frange(find(Pxxrange == max(Pxxrange)));
    Ppeak(aa,:) = max(Pxxrange);
    ps_curves(aa,:) = Pxx;
end

peak_ax = [(Freqpeak(find(Ppeak == max(Ppeak)))) (find(Ppeak == max(Ppeak)))];
Fpeak = round(peak_ax(1));


if (Fpeak-2) >= 1
    [afilt, bfilt] = butter(2, [(Fpeak-2)/(0.5*out.samplerate) (Fpeak+2)/(0.5*out.samplerate)], 'bandpass'); %15
else
    [afilt, bfilt] = butter(2, [(1)/(0.5*out.samplerate) (Fpeak+2)/(0.5*out.samplerate)], 'bandpass'); %15
end

[dd, cc] = butter(2, [0.5/(0.5*out.samplerate) (6)/(0.5*out.samplerate)], 'bandpass'); %15

[b,a]=butter(2,[0.8/(0.5*out.samplerate) ],'low'); %15


for i=1:3
    out.filt{iii,co}(i,:)=filtfilt(afilt,bfilt,data_t(i,:));
    out.l_filt{iii,co}(i,:)=filtfilt(b,a,data_t(i,:));
    out.w_filt{iii,co}(i,:)=filtfilt(dd,cc,data_t(i,:));
end

out.data_t{iii,co}=data_t;

end