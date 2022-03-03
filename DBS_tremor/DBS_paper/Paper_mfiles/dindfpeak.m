function [peak,Pxx,F,cut]=dindfpeak(si,samplerate2)

[a,b]=  butter(3,[0.5/(0.5*samplerate2)], 'high'); %15
filt_cs=(filtfilt(a,b,si));
[Pxx,F] = pwelch(filt_cs, samplerate2, [], floor(samplerate2)*2, samplerate2);

mini=2;
maxi=10;

cut(1)=(find(F==mini));
cut(2)=(find(F==maxi));
dum=Pxx(cut(1):cut(2));
dumf=F(cut(1):cut(2));
peak=dumf(find(dum==max(dum)));

% figure(2)
% plot(F,Pxx)
% hold on
% xlim([0 50])

end