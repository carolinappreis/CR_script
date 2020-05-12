function [s]=zfiltenv(data,afilt,bfilt,co,iii)


for i=1:size(data,1)
    if co==1
        ns_mat={[1 2 3];[1 2 3];[3 2 1];[1 2 3];[3 2 1];[3 2 1];[3 2 1];[3 2 1];[1 2 3];[1 2 3]};
        s.raw(i,:)=data(ns_mat{iii,1}(i),:);
    else
        s.raw(i,:)=data(i,:);
    end
    s.filt(i,:)=filtfilt(bfilt,afilt,s.raw(i,:))*10*9.81/0.5;
%     s.filt(i,:)=filtfilt(bfilt,afilt,s.raw(i,:));
    s.env(i,:)=abs(hilbert(s.filt(i,:)));
    s.phase(i,:)=angle(hilbert(s.filt(i,:)));
    s.freq(i,:)=(smooth((1000/(2*pi))*diff(unwrap(s.phase(i,:))),500))';
    
    s.z(i,:)=zscore(s.raw(i,:));
    s.zfilt(i,:)=filtfilt(bfilt,afilt,s.z(i,:))*10*9.81/0.5;
%     s.zfilt(i,:)=filtfilt(bfilt,afilt,s.z(i,:));
    s.zenv(i,:)=abs(hilbert(s.zfilt(i,:)));
    s.zphase(i,:)=angle(hilbert(s.zfilt(i,:)));
    s.zfreq(i,:)=(smooth((1000/(2*pi))*diff(unwrap(s.zphase(i,:))),500))';
       
end
end


