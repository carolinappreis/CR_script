function [out]=psd_rs(out)
%  load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/ns_seg_pxx.mat')

samplerate=1000;
for iii=1:size(out.start_c,1)
    m_ax=1;
    ARC=squeeze(out.change_c{iii,2}{m_ax,1});
    sup1=[];
    amp1=[];
    
    for y1=1:size(ARC,2)
        for y2=1:size(out.z_seg{iii,y1},1)
            dum=(out.z_seg{iii,y1}(y2,:));
            if ARC(y2,y1)>0
                amp1=[amp1  dum];
            else
                sup1=[sup1  dum];
            end
            clear dum
        end
    end
    
    %     amp=padarray(amp1,[0 2000],0,'both');
    %     sup=padarray(sup1,[0 2000],0,'both');
    
    
    
    [Pxx_s,F]=pwelch(sup1',samplerate*4,[],2*samplerate,samplerate);
    
    [Pxx_a,F]=pwelch(amp1',samplerate*4,[],2*samplerate,samplerate);
    
    [Pxx_ns,F]=pwelch(out.ns(iii,1),samplerate*4,[],2*samplerate,samplerate);
    
    
    
    %     sig=[2 4 5 7 10];
    %     if ismember(iii,sig)
    %         [Pxx_pls,F]=pwelch(out.pls_hu{iii,1},samplerate*4,[],2*samplerate,samplerate);
    %         [Pxx_nsh,F]=pwelch(out.ns_hu{iii,1},samplerate*4,[],2*samplerate,samplerate);
    %
    %     else
    %         Pxx_pls=NaN(1,length(F));
    %         Pxx_nsh=NaN(1,length(F));
    %
    %
    %     end
    
    out.Pxx_s_all(iii,:)=Pxx_s;
    out.Pxx_a_all(iii,:)=Pxx_a;
    out.Pxx_ns_all(iii,:)=Pxx_ns;
   
    %     out.Pxx_pls_all(iii,:)=Pxx_pls;
    %     out.Pxx_nsh_all(iii,:)=Pxx_nsh;
    
    out.F=F;
    %     out.Pxx_ns_all(y,:)=Pxx_ns;
    
    %         clear Pxx_nstim
    %     s{iii,1}=sup1;
    %     a{iii,1}=amp1;
    %     n_s{iii,1}=(ns{iii,1})';
    clear Pxx_a Pxx_s Pxx_ns Pxx_pls Pxx_nsh

    
end
end

