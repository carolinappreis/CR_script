function   [ffr,aam]=simp_new(ax,end3,begin3,freqi,phase,envelope)

% tremor_f2(1,:)=unwrap(phase(ax,begin3:end3));
% tremor_f22(1,:)=(phase(ax,begin3)+(0:1:(end3-begin3))*2*pi/(1000./mean(freqi(ax,begin3-1000:begin3))));
% ffr= (tremor_f2(end)-tremor_f22(end))/(2*pi*0.001*(end3-begin3)); %mean(frequency(end3-1000:end3));%

ffr=NaN;

aam=(mean(envelope(ax,end3-1000:end3))-mean(envelope(ax, begin3-1000:begin3)))./mean(envelope(ax, begin3-1000:begin3));

end



