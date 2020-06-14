function [s]=zfiltenv(d,afilt,bfilt,co,iii,s)
data=d.data_ds;

for i=1:size(data,1)
    if co==1
        ns_mat={[1 2 3];[1 2 3];[3 2 1];[1 2 3];[3 2 1];[3 2 1];[3 2 1];[3 2 1];[1 2 3];[1 2 3]};
        s.raw{iii,co}(i,:)=data(ns_mat{iii,1}(i),:);
    else
        s.raw{iii,co}(i,:)=data(i,:);
    end
    %   s.filt{iii,co}(i,:)=filtfilt(bfilt,afilt,s.raw{iii,co}(i,:))*10*9.81/0.5;
    s.filt{iii,co}(i,:)=filtfilt(bfilt,afilt,s.raw{iii,co}(i,:));
    s.env{iii,co}(i,:)=abs(hilbert(s.filt{iii,co}(i,:)));
    s.phase{iii,co}(i,:)=angle(hilbert(s.filt{iii,co}(i,:)));
    s.freq{iii,co}(i,:)=(smooth((1000/(2*pi))*diff(unwrap(s.phase{iii,co}(i,:))),500))';
    
    s.z{iii,co}(i,:)=zscore(s.raw{iii,co}(i,:));
    %   s.zfilt{iii,co}(i,:)=filtfilt(bfilt,afilt,s.z{iii,co}(i,:))*10*9.81/0.5;
    s.zfilt{iii,co}(i,:)=filtfilt(bfilt,afilt,s.z{iii,co}(i,:));
    s.zenv{iii,co}(i,:)=abs(hilbert(s.zfilt{iii,co}(i,:)));
    s.zphase{iii,co}(i,:)=angle(hilbert(s.zfilt{iii,co}(i,:)));
    s.zfreq{iii,co}(i,:)=(smooth((1000/(2*pi))*diff(unwrap(s.zphase{iii,co}(i,:))),500))';
end

end


