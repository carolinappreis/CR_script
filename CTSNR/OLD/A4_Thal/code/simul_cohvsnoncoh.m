clear all
Fs=1000;
t=0:(1/Fs):(100000-1)*(1/Fs);
y = sin(2*pi*20*t);
x=rand(length(y),1)';
plot(t,y)
hold on
plot(t,x)
filtrange=20;
[b,a]=butter(2,[(filtrange-5)/(0.5*Fs) (filtrange+5)/(0.5*Fs)],'bandpass');
filtrand=filtfilt(b,a,x);
[Pxx_ind,F_ind]=mscohere(x,y,Fs,[],Fs,Fs);
P=sum(Pxx_ind(filtrange-5:filtrange+5))./sum(Pxx_ind(1:end))
dif_angs=angle(hilbert(y))-angle(hilbert(x));

psi=abs(sum(exp(sqrt(-1)*(dif_angs)))./length(dif_angs))



