function   [mm]=simp(ax,iii,end3,begin3,freqi,phase,ns_mat)
tremor_f2(1,:)=unwrap(phase(ns_mat(iii,ax),begin3:end3));
tremor_f22(1,:)=(phase(ns_mat(iii,ax),begin3)+(0:1:(end3-begin3))*2*pi/(1000./mean(freqi(ns_mat(iii,ax),begin3-1000:begin3))));
mm= (tremor_f2(end)-tremor_f22(end))/(2*pi*0.001*(end3-begin3)); %mean(frequency(end3-1000:end3));%
end



