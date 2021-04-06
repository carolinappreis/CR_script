function [pref_in, pref_out, mean_ang, spike_rate]=pref_lock_ang(data,ecog,name)

srn=1000;
[b,a]=butter(2,[15/(0.5*srn) 35/(0.5*srn)],'bandpass');
tt=0;
r=0;
for j=1:size(data,1)
    clearvars -except j b a ecog data srn tt name pref_r srate in_b n_spk r in_bs in_mb in_ms
    epochs_t=[];
    ctx=ecog(j,:);
    %%% ecog_shuf=suffle_data(ctx);   optional surrogate
    Ecogfiltered=filtfilt(b,a,ctx);
    env=abs(hilbert(Ecogfiltered));
    
    [onset,offset]=bursts_aligned(env,ctx);
    onset1=onset; clear onset;  onset=horzcat(onset1{:});
    offset1=offset; clear offset;  offset=horzcat(offset1{:});
    
    data_ones=find(data(j,:)==1);
    hp=wrapToPi(angle(hilbert(Ecogfiltered)));
    ang=hp(data_ones);
    
    data_zeros=find(data(j,:)==0); dat_b=hp;
    dat_b(data_zeros)=NaN;
    
    
    if (circ_rtest(ang))<0.05
        tt=tt+1;
        pref_r(tt,:)=circ_mean(ang');
        
        d=data(j,:);
        srn=1000;
        spkrate_1=[];
        for i =1:srn:(length(d)-srn);
            spkrate_1=[spkrate_1 numel(find(d(i:i+srn)==1))];
        end
        srate(tt,1)=mean(spkrate_1);
        
        for jj=1:length(onset)
            r=r+1;
            el=200;
            if onset(jj)>el && onset(jj)+el<length(dat_b)
                epoch=dat_b(onset(jj)-el:onset(jj)+el);
                dum=find(~isnan(epoch));
                if (numel(dum))>2
                    epochs_t=[epochs_t  epoch(dum)];
                end
            n_spk(1,r)= numel(dum);
            end
        end
        in_b(j,1)=circ_mean(epochs_t');
        in_mb(j,1)=circ_r(epochs_t');
        
        ep_sur=[];
        for n=1:(length(dat_b)/1000)
            idx_sur=randi([el+1,(length(dat_b)-el)],1,1);      
            epoch=dat_b(idx_sur-el:idx_sur+el);
                dum=find(~isnan(epoch));
                if (numel(dum))>2
                    ep_sur=[ep_sur epoch(dum)];
                end
        end
        
        in_bs(j,1)=circ_mean(ep_sur');
        in_ms(j,1)=circ_r(ep_sur');
        
    end
end

pref_mean=rad2deg(circ_mean(pref_r));
if pref_mean<0
    pref_mean=360+pref_mean;
end

% pref_std=rad2deg(circ_std(pref_r));
spike_rate=[mean(srate) std(srate)];

histogram(n_spk)
box('off')
xlabel('number of spikes per bursts')
ylabel('counts all bursts')

pref_in=rad2deg(circ_mean(in_b))
if pref_in<0
    pref_in=360+pref_in;
end
pref_out=rad2deg(circ_mean(in_bs))
if pref_out<0
    pref_out=360+pref_out;
end

mean_ang=[mean(in_mb) mean(in_ms)];
end