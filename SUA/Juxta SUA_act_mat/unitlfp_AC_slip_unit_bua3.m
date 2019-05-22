unitindex=[];
if size(WaveDatalfp_DC,1)>1
    unitindex=2:size(WaveDatalfp_DC,1);
end

if ~isempty(unitindex)
    for i=1:length(unitindex)
        eval(['gpibua' int2str(i) '=WaveDatalfp_DC(' int2str(unitindex(i)) ',:);']);
        eval(['gpiunit' int2str(i) '=WaveData_DC(' int2str(unitindex(i)) ',:);']);
    end
end
clear i

if ~isempty(unitindex)
    for ii=1:length(unitindex)
        eval(['dummy_bua=gpibua' int2str(ii) ';']);
        eval(['dummy_unit=gpiunit' int2str(ii) ';']);
        xr=find(dummy_unit==1);
        if length(xr)>=10
            
            kkk=1;
            
            
            [Pxx_ind,F_ind]=pwelch(dummy_bua,samprate,[],samprate,samprate);
            [Pxx(unitcounter,:),F]=mscohere(WaveDatalfp_DC(1,:),dummy_bua,samprate,[],samprate,samprate);
            Pxx_ind_beta=Pxx(unitcounter,16:36);
            filtrange2=14+find(Pxx_ind_beta==max(Pxx_ind_beta));
            [b,a]=butter(2,[(filtrange2-5)/(0.5*samprate) (filtrange2+5)/(0.5*samprate)],'bandpass');
            clear lfp_or dumdum amp_lfp2
            lfp_ctx=filtfilt(b,a,WaveDatalfp_DC(1,:));
            lfp_or= filtfilt(b,a,dummy_bua);
            
            slope_altered_units_identified
            
            ctx_lfpangle=angle(hilbert(lfp_ctx));
            dumdumctx= hilbert(lfp_ctx);
            amp_ctx=((sqrt((real(dumdumctx)).^2+((imag(dumdumctx)).^2)))');
            ampdumctx=((amp_ctx-(median(amp_ctx)))./(median(amp_ctx))*100)';
            
            gc=ctx_lfpangle;
            ctx_ang_unit(1:length(dummy_unit))=-10;
            ctx_ang_unit(find(dummy_unit==1))=gc(find(dummy_unit==1));
            S=sin(gc);
            
            tral(unitcounter) = circ_rtest(gc(find(dummy_unit==1)));
            vint=NaN(1,length(dummy_unit));
            vasehist(unitcounter,:)=hist(gc(dummy_unit==1),-pi:2*pi/12:pi);
            Ar=find(dummy_unit==1)./samprate;
            vr=1./diff(Ar);
            
            vint = interp1(xr(2:end),vr,xr(2):1:length(gpiunit1),'next');
            vint = [NaN(1,xr(2)-1) vint];
            freq_unit(unitcounter)=nanmedian(vint);%(1000/(2*windowsize))*
            
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
            freq_ac_unit=NaN(l_duration,80,100);
            for i=2:l_duration-1
                
                if beginindex22(i)-2050>0 && (endindex22(i)+2050)<length(gc) ...
                        && beginindex22(i)-endindex22(i-1)>200
                    beginindex3=randi([2051 length(gc)-2051],1,1);
                    iii=-2000:50:1950;
                    for iii_index=1:80
                        
                        ind_a=round(abs(gc(beginindex22(i)+iii(iii_index)-30:beginindex22(i)+iii(iii_index)+30))*10)/10;
                        ind_b=round(abs(gc(beginindex22(i)+iii(iii_index)+20:beginindex22(i)+iii(iii_index)+80))*10)/10;
                        a=max(find(ind_a==min(ind_a)))+beginindex22(i)+iii(iii_index)-31;
                        b=max(find(ind_b==min(ind_b)))+beginindex22(i)+iii(iii_index)+19;
                        update=10;
                        while a==b
                            ind_b=round(abs(gc((beginindex22(i)+iii(iii_index)+update+20):(beginindex22(i)+iii(iii_index)+80+update)))*10)/10;
                            a=max(find(ind_a==min(ind_a)))+beginindex22(i)+iii(iii_index)-31;
                            b=max(find(ind_b==min(ind_b)))+beginindex22(i)+iii(iii_index)+19+update;
                            update=update+10; 
                        end
                        
                        update=10;
                        while b-a<25
                            ind_b=round(abs(gc((beginindex22(i)+iii(iii_index)+update+20):(beginindex22(i)+iii(iii_index)+80+update)))*10)/10;
                            a=max(find(ind_a==min(ind_a)))+beginindex22(i)+iii(iii_index)-31;
                            b=max(find(ind_b==min(ind_b)))+beginindex22(i)+iii(iii_index)+19+update;
                            update=update+10;   
                        end
                        
                        update=10;
                        while b-a>70
                            ind_b=round(abs(gc((beginindex22(i)+iii(iii_index)-update+20):(beginindex22(i)+iii(iii_index)+80-update)))*10)/10;
                            a=max(find(ind_a==min(ind_a)))+beginindex22(i)+iii(iii_index)-31;
                            b=max(find(ind_b==min(ind_b)))+beginindex22(i)+iii(iii_index)+19-update;
                            update=update+10;  
                        end

                        
                        %if ~isempty(find(amp_lfp3(beginindex22(i):endindex22(i))==1))
                        freq_ac(kkk,iii_index)=sum(slips3(a:b));
                        %freq_ac5(kkk,iii,:)=phase_dif(a:a+50);
                        freq_ac_unit(kkk,iii_index,1:(b-a+1))=dummy_unit(a:b);%ctx_ang_unit(a:b);
                        %end
                        clear ind_a a X ind_b b
                    end
                    dur_ac(kkk)=duration22(i);
                    kkk=kkk+1;
                    clear X Y
                end
            end
            clear i
            
            if unitcounter==33
                hc=1;
            end
            %
            %             k=1;
            %             for i=1:50:3951
            %                 freq_ac55(:,k)=sum(freq_ac(:,i:i+50));
            %                 k=k+1;
            %             end
            
            freq_ac2(ii,1:80)= nanmean(freq_ac,1)/(size(freq_ac,1));
            
            for i=1:80
                dumvar1(1:size(freq_ac_unit,1),1:size(freq_ac_unit,3))=freq_ac_unit(:,i,:);
                dumvar2=dumvar1(:);
                if length(find(dumvar2~=-10))>=(size(freq_ac,1)/10) %~isempty(find(dumvar2~=-10))
                    dumvar3=dumvar2(find(dumvar2~=-10 & ~isnan(dumvar2)));
                    freq_ac2_unit(ii,i)=abs(nansum(exp(sqrt(-1).*dumvar3))./(length(dumvar3)));%circ_rtest(dumvar3);
                else
                    freq_ac2_unit(ii,i)=NaN;
                end
                clear dumvar1 dumvar2 dumvar3
            end
            
        else
            freq_ac2(ii,1:80)=NaN;
            freq_ac2_unit(ii,1:80)=NaN;
            tral(unitcounter) = 1;
            freq_unit(unitcounter) = NaN;
        end
        unitcounter=unitcounter+1;
    end
end

