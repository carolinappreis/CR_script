clear all
% d=load('/Users/Carolina/Documents/NIC/20211118175909_Patient01_pilot_no brain.easy');
d=load('/Users/Carolina/Documents/NIC/20211118175252_Patient01_pilot_no brain.easy');
% d=load('/Users/Carolina/Documents/NIC/20130925192801_PatientTMS.easy');

% NE_Viewer('20130925192801_PatientTMS.easy',999,500,{'Fp1','Fp2','F3','F4','C3','C4','P3','Cz'});

% NE_Viewer('20211118175252_Patient01_pilot_no brain.easy',999,500,{'Fpz','Fz','C3','Cz','C4','O1','Oz','O2'});

samplerate=500;
cond_t={(1:2*60*samplerate);(2*60*samplerate:3.5*60*samplerate)};
cond_label={'Rest EO';'Rest EC'};

for cc=1:2
    for ii=1:8
        % [Pxx,F] = pwelch(d((2*60*500):(3.5*60*500),ii),samplerate, [], samplerate, samplerate);
        [Pxx,F] = pwelch(d((cond_t{cc,1}),ii),2*samplerate, [], 2*samplerate, samplerate);        
        power(cc,ii,:)=Pxx; clear Pxx        
    end
end


% subplot(8,1,ii)
% xlim(log10([4 48]))
plot(log10(F),log10(squeeze(power(2,7,:))))
set(gca,'XScale','log','YScale','log')
new=squeeze(power(1,:,:));

plot(F,squeeze(new'))
xlim([0 50])

[p,q]=butter(3,[(0.5/(0.5*samplerate))],'high');
fi=filtfilt(p,q,d(:,6));
