
function [sua] = makemua_CR_spikes(hpsig,msback,msfoward,muasr,rssr,uthresh)

% hpsig = wideband signal
% msback = length of data replaced before a spike in s
% msback = length of data replaced before after spike in s
% muasr = sampling rate of the orignal channel
% rssr = sampling rate of the new channel
% uthresh = threshold for single units in SDs of baseline. 


% For example
%[mua] = makemua_fede(muatemp,0.001,0.003,muasr,4);
% 
% hpsig = unit (already high pass filtered);
% msback = 0.001;
% msfoward = 0.003;
% muasr = 16000;
% uthresh = 4;
% Highpass filter signal

clear mua
clear beta*
clear muatemphp
[betab,betaa] = butter(3,300/(muasr/2),'high');
muatemphp = filtfilt(betab,betaa,hpsig);


% clear muatemphp
% muatemphp = hpsig;

clear E1d*
clear E1n*
% Resample 
E1d = muatemphp;
tback = floor(msback*muasr);
tforward = floor(msfoward*muasr);
clear E1 abs;
E1dabs = abs(E1d);

clear tthres*
% Set threshold for single units
tthres = mean(E1dabs)+(std(E1dabs)*uthresh);
% tthresneg = mean(E1d)-(std(E1d)*uthresh);

clear tpass*

% Find index of units
tpass = find(E1dabs>tthres);
% tpassneg = find(E1d<tthresneg);

% Make binary unit channel for pos and neg crosses 
tpassbin = zeros(size(E1d,1),1);
tpassbin(tpass) = 1;


srn=rssr;
dataold=tpassbin;
dataold=full(dataold);
data=zeros(1,100000);
timeold=0:1/muasr:(size(dataold,1)-1)/muasr;
time=0:1/srn:(size(data,2)-1)/srn;
spk_t=timeold(tpass);
spk_tround=round(spk_t,3);
n=[];
for i=1:length(spk_t)
    [ d, ix ] = min( abs( time-spk_tround(i) ) );
    n=[n ix];
end
data(n)=1;

clear sua 
sua = data;
end



% 
% clear mua
% clear tx1
% mua =E1z;


% [pow,f] = pwelch(mua,1000,500,1000,1000);



% figure(1)%mix of noise (all noise, instead of data PLUS noise)
% plot(E1z,'r')
% xlim([0 500])
% ylim([-0.05 0.05])
% figure(2)%noise plus data
% plot(E1d)
% xlim([0 500])
% ylim([-0.05 0.05])
% figure(3)
% plot(allnoisematch)
% xlim([0 500])
% ylim([-0.05 0.05])
% figure(4)% segments of data are seen in blue among noise (red)
% plot(E1d)
% hold on
% plot(E1z,'r')
% hold on
% xlim([0 500])
% ylim([-0.05 0.05])