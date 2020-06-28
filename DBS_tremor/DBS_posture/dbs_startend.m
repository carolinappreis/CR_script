function [peak_ax, start, ending, yy, h_up]= dbs_startend(d,samplerateold,samplerate,iii,co,start,ending,yy,h_up)

data=d.data_raw;  tremor_ds=d.data_ds;

if co==1
    
    NS_BE_P
    start{iii,co}=hu{iii,:};
    ending{iii,co}=hd{iii,:};
    
    
    % when patient's hand is up
    handup = [];
    for i = 1:length(start{iii,co})
        handup = [handup start{iii,co}(i):ending{iii,co}(i)]; %#ok<*AGROW>
    end
    clear i
    handup = sort(handup,'ascend');
    
    h_up{iii,co}=handup;
    
    for aa = 1:3
        [Pxx,F] = pwelch(tremor_ds(aa,handup), samplerate, [], samplerate, samplerate);
        %          [Pxx,F] = pwelch(tremor_ds(aa,:), samplerate, [], samplerate, samplerate);
        
        frange = F(3:10);
        Pxxrange = Pxx(3:10);
        Freqpeak(aa,:) = frange(find(Pxxrange == max(Pxxrange)));
        Ppeak(aa,:) = max(Pxxrange);
        ps_curves(aa,:) = Pxx;
    end
    
    peak_ax = [(Freqpeak(find(Ppeak == max(Ppeak)))) (find(Ppeak == max(Ppeak)))];
    Fpeak = peak_ax(1);
    
    
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
    indexes4 = [index(1) index(find(diff(index) ./ samplerateold > 0.95)+1)];
    indexes3 = [index(find(diff(index) ./ samplerateold > 0.95)) index(end)];
    
    
    dd2 = round(data(4, :) *100) ./ 100;
    for i = 1:length(indexes4)
        xx(i) = round(dd2(indexes4(i)) ./ 0.1); %#ok<*SAGROW>
    end
    clear i
    
    
    ss = floor((indexes4 ./ samplerateold)*samplerate)+addon;
    ee = floor((indexes3 ./ samplerateold)*samplerate)+addon+addon_end;%floor(5*samplerate);
    
    % when patient's hand is up
    handup = [];
    for i = 1:length(ss)
        handup = [handup ss(i):ee(i)]; %#ok<*AGROW>
    end
    clear i
    handup = sort(handup,'ascend');
    h_up{iii,co}=handup;
    
    for aa = 1:3
        [Pxx,F] = pwelch(tremor_ds(aa,handup), samplerate, [], samplerate, samplerate);
        %          [Pxx,F] = pwelch(tremor_ds(aa,:), samplerate, [], samplerate, samplerate);
        
        frange = F(3:10);
        Pxxrange = Pxx(3:10);
        Freqpeak(aa,:) = frange(find(Pxxrange == max(Pxxrange)));
        Ppeak(aa,:) = max(Pxxrange);
        ps_curves(aa,:) = Pxx;
    end
    
    peak_ax = [(Freqpeak(find(Ppeak == max(Ppeak)))) (find(Ppeak == max(Ppeak)))];
    Fpeak=peak_ax(1);
    
    new = find(data(2,:) > 4);
    difp = find((diff(new)) > 100000); % are you trying to threshold at 9.6 seconds?
    ep_1 = [new(difp) new(end)];
    sp_1 = [new(1) new(difp+1)];
    
    %%% input start all trial
    start_t=1;
    sp=sp_1(1,start_t:end);
    ep=ep_1(1,start_t:end);
   
    
%     time=1:length(data(1,:));
%     plot(time,data(4,:))
%     hold on
%     plot(time(sp),data(4,sp),'r.')
%     plot(time(ep),data(4,ep),'k.')
    
    if co==2
        %%% posture/spiral trials
        dt1_s=[sp_1(start_t:2:end)];dt1_e=[ep_1(start_t:2:end)];
        dt2_s=sp_1(start_t+1:2:end); dt2_e=ep_1(start_t+1:2:end);
        
        
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
        
        
        %     figure()
        %     time=1:length(data(4,:));
        %     plot(time,data(4,:))
        %     hold on
        %     plot(time(index),data(4,index),'r.')
        %     plot(time(start2),data(4,start2),'ko')
        %     plot(time(ending2),data(4,ending2),'bo')
        
        
        clear start ending xx
        pstart{1,1}=floor((start1./samplerateold)*samplerate)+addon;
        pending{1,1}=floor((ending1./samplerateold)*samplerate)+addon+addon_end;%floor(5*samplerate);
        pstart{2,1}=floor((start2./samplerateold)*samplerate)+addon;
        pending{2,1}=floor((ending2./samplerateold)*samplerate)+addon+addon_end;%floor(5*samplerate);
        
        xx{1,1}=xx1;
        xx{2,1}=xx2;
        
        %%% choosing start{1,1}/ending{1,1}/xx{1,1} to get posture only _ check!
        
        start{iii,co}= pstart{1,1};
        ending{iii,co} = pending{1,1};
        yy{iii,co} = xx{1,1};
        
        
        %         figure()
        %         index_ds=floor((index./samplerateold)*samplerate)+addon;
        %         time=1:length(tremor_ds);
        %         plot(time,tremor_ds(1,:))
        %         hold on
        %         plot(time(1,index_ds),tremor_ds(1,index_ds),'r.')
        %         plot(time(start{iii,co}),tremor_ds(1,start{iii,co}),'bo')
        %         plot(time(ending{iii,co}),tremor_ds(1,ending{iii,co}),'ko')
        %----------------------------------------------------------------------------------------------
    elseif co==3 
        start{iii,co}= floor((sp./samplerateold)*samplerate)+addon;
        ending{iii,co} = floor((ep./samplerateold)*samplerate)+addon+addon_end;
    end
    
end
end


