clear all
d=load('/Users/Carolina/Documents/NIC/20211118175909_Patient02_pilot_no brain.easy');
% d=load('/Users/Carolina/Documents/NIC/20130925192801_PatientTMS.easy');

% NE_Viewer('20130925192801_PatientTMS.easy',999,500,{'Fp1','Fp2','F3','F4','C3','C4','P3','Cz'});

% NE_Viewer('20211118161511_Patient01.easy',999,500,{'Fp1','Fpz','Fp2','Fz','C3','Cz','C4','Pz'});

samplerate=500;


for ii=1:8
% [Pxx,F] = pwelch(d((2*60*500):(3.5*60*500),ii),samplerate, [], samplerate, samplerate);
[Pxx,F] = pwelch(d((3.4*60*500):end),ii),2*samplerate, [], 2*samplerate, samplerate);

power(ii,:)=Pxx; clear Pxx
% subplot(8,1,ii)


% xlim(log10([4 48]))
end


plot(log10(F),log10(power(5,:)))
set(gca,'XScale','log','YScale','log')

plot(F,power(5,:))
xlim([0 50])