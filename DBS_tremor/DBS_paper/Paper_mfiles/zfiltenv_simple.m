function [s]=zfiltenv_simple(d,peak_ax,co,iii,s,samplerate)

data=d.data_ds;
Fpeak=peak_ax(1);

if (Fpeak-2) >= 1
    [afilt, bfilt] = butter(2, [(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)], 'bandpass'); %15
else
    [afilt, bfilt] = butter(2, [(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)], 'bandpass'); %15
end

    [h1, h2] = butter(3,[1/(0.5*d.samplerateold)],'high');


for i=1:size(data,1)
    s.raw{iii,co}(i,:)=data(i,:);
    s.filt{iii,co}(i,:)=filtfilt(afilt,bfilt,s.raw{iii,co}(i,:));
    s.filth{iii,co}(i,:)=filtfilt(h1,h2,s.raw{iii,co}(i,:));
    s.env{iii,co}(i,:)=abs(hilbert(s.filt{iii,co}(i,:)));
    s.env_acc{iii,co}(i,:)=abs(hilbert(s.filt{iii,co}(i,:).*(10*9.81/0.5)));
    s.phase{iii,co}(i,:)=angle(hilbert(s.filt{iii,co}(i,:)));
    s.freq{iii,co}(i,:)=(smooth((1000/(2*pi))*diff(unwrap(s.phase{iii,co}(i,:))),500))';
    
    s.z{iii,co}(i,:)=zscore(data(i,:));
    s.zfilt{iii,co}(i,:)=filtfilt(afilt,bfilt,s.z{iii,co}(i,:));
    s.zfilth{iii,co}(i,:)=filtfilt(h1,h2,s.z{iii,co}(i,:));
    s.zenv{iii,co}(i,:)=abs(hilbert(s.zfilt{iii,co}(i,:)));
    s.zphase{iii,co}(i,:)=angle(hilbert(s.zfilt{iii,co}(i,:)));
    s.zfreq{iii,co}(i,:)=(smooth((1000/(2*pi))*diff(unwrap(s.zphase{iii,co}(i,:))),500))';
end

%%% for combined signal
% data_2=sqrt(data(1,:).^2+data(2,:).^2+data(3,:).^2);
% data_3=repmat(sqrt(data(1,:).^2+data(2,:).^2+data(3,:).^2),3,1);
% 
% 
% for i=1:size(data_3,1)
%     s.raw{iii,co}(i,:)=data_3(i,:);
%     s.filt{iii,co}(i,:)=filtfilt(afilt,bfilt,s.raw{iii,co}(i,:));
%     s.filth{iii,co}(i,:)=filtfilt(h1,h2,s.raw{iii,co}(i,:));
%     s.env{iii,co}(i,:)=abs(hilbert(s.filt{iii,co}(i,:)));
%     s.env_acc{iii,co}(i,:)=abs(hilbert(s.filt{iii,co}(i,:).*(10*9.81/0.5)));
%     s.phase{iii,co}(i,:)=angle(hilbert(s.filt{iii,co}(i,:)));
%     s.freq{iii,co}(i,:)=(smooth((1000/(2*pi))*diff(unwrap(s.phase{iii,co}(i,:))),500))';
%     
%     s.z{iii,co}(i,:)=zscore(data_3(i,:));
%     s.zfilt{iii,co}(i,:)=filtfilt(afilt,bfilt,s.z{iii,co}(i,:));
%     s.zfilth{iii,co}(i,:)=filtfilt(h1,h2,s.z{iii,co}(i,:));
%     s.zenv{iii,co}(i,:)=abs(hilbert(s.zfilt{iii,co}(i,:)));
%     s.zphase{iii,co}(i,:)=angle(hilbert(s.zfilt{iii,co}(i,:)));
%     s.zfreq{iii,co}(i,:)=(smooth((1000/(2*pi))*diff(unwrap(s.zphase{iii,co}(i,:))),500))';
% 
% end


end


