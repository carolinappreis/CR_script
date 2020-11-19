unitindex=[];
if size(WaveData_DC,1)>1
    unitindex=2:size(WaveData_DC,1);
end

if ~isempty(unitindex)
    for i=1:length(unitindex)
        eval(['gpiunit' int2str(i) '=WaveData_DC(' int2str(unitindex(i)) ',:);']);
    end
end

if ~isempty(unitindex)
    for ii=1:length(unitindex)
        kkk=1;
        eval(['dummy=gpiunit' int2str(ii) ';']);
        
        [Pxx_ind,F_ind]=pwelch(dummy,samprate,[],samprate,samprate);
        [Pxx(unitcounter,:),F]=mscohere(WaveData_DC(1,:),dummy,samprate,[],samprate,samprate);
        Pxx_ind_beta=Pxx(unitcounter,16:36);
        filtrange2=14+find(Pxx_ind_beta==max(Pxx_ind_beta));
        [b,a]=butter(2,[(filtrange2-5)/(0.5*samprate) (filtrange2+5)/(0.5*samprate)],'bandpass');
        clear lfp_or dumdum amp_lfp2
        lfp_ctx=filtfilt(b,a,WaveData_DC(1,:));
        lfp_or= filtfilt(b,a,dummy);
        
        slope_altered_units_identified
        
        ctx_lfpangle=angle(hilbert(lfp_ctx));
        dumdumctx= hilbert(lfp_ctx);
        amp_ctx=((sqrt((real(dumdumctx)).^2+((imag(dumdumctx)).^2)))');
        ampdumctx=((amp_ctx-(median(amp_ctx)))./(median(amp_ctx))*100)';
        
        %         ctx_lfpangle=angle(hilbert(lfp_ctx));
        %         dumdumctx= hilbert(lfp_ctx);
        %         amp_ctx=((sqrt((real(dumdumctx)).^2+((imag(dumdumctx)).^2)))');
        %         [b,a]=butter(2,[20/(0.5*samprate) 30/(0.5*samprate)],'bandpass');
        %         noise_ctx=filtfilt(b,a,randn(size(time)));
        %         noise_angle=angle(hilbert(noise_ctx));
        %         lfp_ctx=amp_ctx'.*cos(noise_angle);
        %         ctx_lfpangle=angle(hilbert(lfp_ctx));
        %         dumdumctx= hilbert(lfp_ctx);
        %         ampdumctx=((amp_ctx-(median(amp_ctx)))./(median(amp_ctx))*100)';
        
        gc=ctx_lfpangle;%angle(dumdum)-angle(dumdumstn); %
        S=sin(gc);
        
        duration22=duration2;%(find(duration2>median(duration2))); %duration>= 50 & duration<= 200
        beginindex22=beginindex2;%(find(duration2>median(duration2)));
        endindex22=endindex2;%(find(duration2>median(duration2)));
        l_duration=length(duration22);
        time=0:0.001:(length(gc)-1)/1000;
        
        
        
        phase_1=angle(hilbert(lfp_or));
        env=abs(hilbert(lfp_or));
        env2=((env-(median(env)))./(median(env))*100)';
        phase_2=angle(dumdumctx);
        instphase=diff(unwrap(phase_1));
        freq=medfilt1((1000/(2*pi))*instphase,51);
        instphase2=diff(unwrap(phase_2));
        freq2=medfilt1((1000/(2*pi))*instphase2,51);
        freq3=(freq2-freq);
        freq4=(freq3-nanmean(freq3))./nanstd(freq3);
        phase_dif=phase_1-phase_2;
        phase_dif(phase_dif<-pi)=phase_dif(phase_dif<-pi)+2*pi;
        phase_dif(phase_dif>pi)=phase_dif(phase_dif>pi)-2*pi;
        slips=diff(unwrap(phase_dif));
        slips2=(slips-mean(slips))./std(slips);
        slips3=zeros(1,size(slips2,2));
        slips3(slips2>1.96)=1;
        slips3(slips2<-1.96)=1;
        
        ind_thr=prctile(env,75);
        amp_lfp3=zeros(size(env));
        amp_lfp3(env>ind_thr)=1;
        
        vasehist(unitcounter,:)=hist(phase_dif,-pi:2*pi/12:pi);
        tral(unitcounter) = circ_rtest(phase_dif);
        for i=2:l_duration-1
            iii=0;
            if beginindex22(i)-550>0 && (endindex22(i)+550)<length(gc) ...
                    && beginindex22(i)-endindex22(i-1)>200
                
                beginindex3=randi([2050 length(amp_ctx)-2050],1,1);
                
                ind_a=round(abs(gc(beginindex22(i)+iii-30:beginindex22(i)+iii+30))*10)/10;
                ind_b=(beginindex22(i)-500):(beginindex22(i)+500);
                a=max(find(ind_a==min(ind_a)))+beginindex22(i)+iii-31;
                amp_lfp2(i,:)= (env(a-500:a+500)-median(env(a-500:a)))./median(env(a-500:a));
%                 if isempty(find(amp_lfp3(beginindex22(i):endindex22(i))==1))
                    X=slips3(beginindex22(i)-500:beginindex22(i)+500);
                    XX=phase_dif(beginindex22(i)-500:beginindex22(i)+500);
                    %(freq(beginindex22(i)-500:beginindex22(i)+500))-nanmean(freq(beginindex22(i)-500:beginindex22(i)));
                    %(dummy(a-500:a+500)); %phase_dif(beginindex22(i)-500:beginindex22(i)+500);
                    freq_ac(kkk,:)=X;%((X)-mean(dummy))./std(dummy);%X;%(medfilt1(diff(unwrap(X)),51)); %
                    freq_ac5(kkk,:)=phase_dif(beginindex22(i)-500:beginindex22(i)+500);
                    dur_ac(kkk)=duration22(i);
                    kkk=kkk+1;
%                 end
                clear ind_a a X ind_b b
            end
            clear X Y
        end
        k=1;
        for i=1:10:991
            freq_ac55(:,k)=sum(freq_ac(:,i:i+10));
            k=k+1;
        end
        
        for i=1:size(freq_ac,1);
            if dur_ac(i)<500
                if ~isempty(find(freq_ac(i,101:500)==1))
                    endpoint(i)=max(find(freq_ac(i,101:500)==1))+100; %#ok<*SAGROW>
                    avgangle(i)=angle(nansum(exp(sqrt(-1).*freq_ac5(i,endpoint(i)-50:endpoint(i))))./51)-angle(nansum(exp(sqrt(-1).*freq_ac5(i,500:550)))./51);
                else
                    endpoint(i)=NaN;
                    avgangle(i)=NaN;
                end
            else
                if ~isempty(find(freq_ac(i,101:1000)==1))
                    endpoint(i)=max(find(freq_ac(i,101:949)==1))+100; %#ok<*SAGROW>
                    avgangle(i)=angle(nansum(exp(sqrt(-1).*freq_ac5(i,endpoint(i)-50:endpoint(i))))./51)-angle(nansum(exp(sqrt(-1).*freq_ac5(i,500:550)))./51);
                else
                    endpoint(i)=NaN;
                    avgangle(i)=NaN;
                end
            end
        end
        avgangle(avgangle<-pi)=avgangle(avgangle<-pi)+2*pi;
        avgangle(avgangle>pi)=avgangle(avgangle>pi)-2*pi;
        
        [a1,b1]=sort(dur_ac,'descend');
        if length(b1)>25
            freq_ac2(ii,1:25,1:1001)= freq_ac(b1(1:25),1:1001);%nanmean(freq_ac,1);
        else
            freq_ac2(ii,1:25,1:1001)= NaN(25,1001);
        end
        output33_angle(ii)=angle(nansum(exp(sqrt(-1).*avgangle))./size(avgangle,2));
        output33_psi(ii)=abs(nansum(exp(sqrt(-1).*avgangle))./size(avgangle,2));
        %[angle(nansum(exp(sqrt(-1).*avgangle(endpoint<300))))./sum(~isnan(avgangle(endpoint<300))) angle(nansum(exp(sqrt(-1).*avgangle(endpoint>300))))./sum(~isnan(avgangle(endpoint>300)))] ; 
        %nanmean(freq_ac,1); %abs(nansum(exp(sqrt(-1).*freq_ac))./(k-1));
        %nansum(freq_ac,1)./(k-1);% angle(nansum(exp(sqrt(-1).*freq_ac))./(k-1));
        unitcounter=unitcounter+1;
    end
end

