function [afilt, bfilt, start, ending, yy, h_up]= startend(d,samplerateold,samplerate,iii,co,start,ending,yy,h_up)

data=d.data_raw;  tremor_ds=d.data_ds;
pls_phastim={[8 8 10 10];[5 5 8 8 8 8 8];repmat(8,1,8);[5 5]; repmat(8,1,5);repmat(5,1,5); repmat(2,1,6);[9 9 12 12];9;[2 2 12 12]}; % stim phase-phaselocked
    
if co==1

before_ns

start{iii,co}=round(hu{iii,:});
ending{iii,co}=round(hd{iii,:});


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
    
%     plot(ps_curves')
%     xlim([3 30])
%     box('off')
%     ylabel('Power spectra')
%     xlabel('Frequency (Hz)')
end

peak_ax = [(Freqpeak(find(Ppeak == max(Ppeak)))) (find(Ppeak == max(Ppeak)))];
Fpeak = peak_ax(1);

if (Fpeak-2) >= 1
    [afilt, bfilt] = butter(2, [(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)], 'bandpass'); %15
else
    [afilt, bfilt] = butter(2, [(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)], 'bandpass'); %15
end

yy{iii,co}=NaN;

%----------------------------------------------------------------------------------------------
else
    
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

     
    ss = floor((indexes4 ./ samplerateold)*samplerate);%+addon;
    ee = floor((indexes3 ./ samplerateold)*samplerate);%+addon+addon_end;%floor(5*samplerate);
    
    % when patient's hand is up
    handup = [];
    for i = 1:length(ss)
        handup = [handup ss(i):ee(i)]; %#ok<*AGROW>
    end
    clear i
    handup = sort(handup,'ascend');
    h_up{iii,co}=handup;
    
    [Pxx, F] = pwelch(tremor_ds(1,handup), samplerate, [], samplerate, samplerate);
    
    frange = F(3:10);
    Pxxrange = Pxx(3:10);
    
    Fpeak = frange(find(Pxxrange == max(Pxxrange))); %#ok<*FNDSB>
    
    if (Fpeak-2) >= 1
        [afilt, bfilt] = butter(2, [(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)], 'bandpass'); %15
    else
        [afilt, bfilt] = butter(2, [(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)], 'bandpass'); %15
    end
    
    
    new = find(data(2,:) > 4);
    difp = find((diff(new)) > 100000); % are you trying to threshold at 9.6 seconds?
    ep_1 = [new(difp) new(end)];
    sp_1 = [new(1) new(difp+1)];
    
    %%% input start all trial
    start_t = 1;
    sp = sp_1(1, start_t:end);
    ep = ep_1(1, start_t:end);
    
    for ik = 1:length(sp) %%find double start and end points in a stimulation run
        s = (find(([indexes4 >= sp(ik)] + [indexes4 <= ep(ik)]) == 2));
        e = (find(([indexes3 >= sp(ik)] + [indexes3 <= ep(ik)]) == 2));
        tks = (find(diff(xx(s)) == 0)) + 1;
        tke = (find(diff(xx(e)) == 0));
        
        indexes4(s(tks)) = NaN;
        indexes3(e(tke)) = NaN;
        xx(e(tke)) = NaN;
    end
    
    
    if co==2
        %%%% find runs with trigering issues (too few, too many pulses) %5 sec
        %%%% stim
        th1 = (Fpeak * 5) ./ 2;
        th2 = (Fpeak * 5) + round((Fpeak * 5) ./ 2);
        
    elseif co==3
        %%% find runs with trigering issues (too few, too many pulses) 90 sec
        %%% stim
        th1=(Fpeak*90)./2;
        th2=(Fpeak*90)+round((Fpeak*90)./2);
    end
    
    for it = 1:length(indexes4)
        if numel(index(find(index == indexes4(it)):find(index == indexes3(it)))) >= th1 && numel(index(find(index == indexes4(it)):find(index == indexes3(it)))) <= th2
            indexes4(it) = indexes4(it);
            indexes3(it) = indexes3(it);
            xx(it) = xx(it);
        else
            indexes4(it) = NaN;
            indexes3(it) = NaN;
            xx(it) = NaN;
        end
    end
    
    %%%%%%%%%%%%%%%
    indexes4 = indexes4(~isnan(indexes4));
    indexes3 = indexes3(~isnan(indexes3));
    xx = xx(~isnan(xx));
    
    start1 = [];
    ending1 = [];
    xx1 = [];
    for il = 1:length(sp)
        start1 = [start1 indexes4(find(([indexes4 >= sp(il)] + [indexes4 <= ep(il)]) == 2))]; % intersect([1 2 3],[3 4 5])
        ending1 = [ending1 indexes3(find(([indexes3 >= sp(il)] + [indexes3 <= ep(il)]) == 2))];
        xx1 = [xx1 xx(find(([indexes4 >= sp(il)] + [indexes4 <= ep(il)]) == 2))];
    end
    
    
    start{iii,co}= floor((start1 ./ samplerateold) * samplerate);%+addon;
    ending{iii,co} = floor((ending1 ./ samplerateold) * samplerate);%+addon+addon_end;%floor(5*samplerate);
    clear xx
    

    yy{iii,co} = xx1;

    
    %----------------------------------------------------------------------------------------------
    
end
end

