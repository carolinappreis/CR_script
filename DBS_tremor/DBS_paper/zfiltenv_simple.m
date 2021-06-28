function [s]=zfiltenv_simple(d,peak_ax,co,iii,s,samplerate)

data=d.data_ds;
Fpeak=peak_ax(1);

if (Fpeak-2) >= 1
    [afilt, bfilt] = butter(2, [(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)], 'bandpass'); %15
else
    [afilt, bfilt] = butter(2, [(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)], 'bandpass'); %15
end


for i=1:size(data,1)
    s.raw{iii,co}(i,:)=data(i,:);
    s.filt{iii,co}(i,:)=filtfilt(afilt,bfilt,s.raw{iii,co}(i,:));
    s.env{iii,co}(i,:)=abs(hilbert(s.filt{iii,co}(i,:)));
    s.env_acc{iii,co}(i,:)=abs(hilbert(s.filt{iii,co}(i,:).*(10*9.81/0.5)));
    s.phase{iii,co}(i,:)=angle(hilbert(s.filt{iii,co}(i,:)));
    s.freq{iii,co}(i,:)=(smooth((1000/(2*pi))*diff(unwrap(s.phase{iii,co}(i,:))),500))';
    
    s.z{iii,co}(i,:)=zscore(data(i,:));
    s.zfilt{iii,co}(i,:)=filtfilt(afilt,bfilt,s.z{iii,co}(i,:));
    s.zenv{iii,co}(i,:)=abs(hilbert(s.zfilt{iii,co}(i,:)));
    s.zphase{iii,co}(i,:)=angle(hilbert(s.zfilt{iii,co}(i,:)));
    s.zfreq{iii,co}(i,:)=(smooth((1000/(2*pi))*diff(unwrap(s.zphase{iii,co}(i,:))),500))';
end


end


