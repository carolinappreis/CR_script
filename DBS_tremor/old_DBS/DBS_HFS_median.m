% % clear all
% % iii=[1];
% % 
% % for numb=1;
% %     %     :length(iii);
% %     %     clearvars -except nostimout iii numb cc cond nostim
% %     
% %     
% %     for   in2=1:3;
% %         DBS_Fpeak
% %         %     load(strcat('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\0',num2str(iii(numb)),cond{cc,1}))
% %         load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/0',num2str(iii(numb)),'_HFS_PS.mat'));
% %         segmentb=[1621 112001 205601 287501 366701 450401];
% %         segmente=[76121 179301 263501 350901 422701 517501];
% %         start{1,1}=segmentb(1:2:end);
% %         start{2,1}=segmentb(2:2:end);
% %         ending{1,1}=segmente(1:2:end);
% %         ending{2,1}=segmente(2:2:end);
% %         % analysing the "main tremor axis"
% %         
% %         if in2==1
% %             in=3;
% %         elseif in2==2 % other axis 1
% %             in=6;
% %         elseif in2==3 % other axis 2
% %             in=7;
% %         end
% %         data=SmrData.WvData;
% %         samplerateold=SmrData.SR;
% %         tremor=(data(in,:));
% %         addon=92; addon_end=35;
% %         
% %         time=0:1/samplerateold:(size(data,2)-1)/samplerateold;
% %         
% %         %%% downsample
% %         
% %         ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
% %         ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
% %         tremor2(1:size(ts1.data,3))=ts1.data;
% %         samplerate=1000;
% %         
% %         if (Fpeak-2)>=1
% %             [b,a]=butter(2,[(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
% %         else
% %             [b,a]=butter(2,[(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
% %         end
% %         tremor_or=filtfilt(b,a,tremor2)*10*9.81/0.5;
% %         [b,a]=butter(2,[0.8/(0.5*samplerate) ],'low'); %15
% %         % tremor_or=zscore(tremor_or);
% %         dummy=hilbert(tremor_or);
% %         envelope=sqrt((real(dummy).^2)+(imag(dummy).^2));
% %         phase=angle(dummy);
% %         frequency=(smooth((1000/(2*pi))*diff(unwrap(angle(dummy))),500))';
% %         tremor=(data(3,:));% %score(:,1)';%
% %         ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
% %         ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
% %         tremorx(1:size(ts1.data,3))=ts1.data;
% %         filt_x=filtfilt(b,a,tremorx);
% %         tremor=(data(6,:));% %score(:,1)';%
% %         ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
% %         ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
% %         tremory(1:size(ts1.data,3))=ts1.data;
% %         filt_y=filtfilt(b,a,tremory);
% %         tremor=(data(7,:));% %score(:,1)';%
% %         ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
% %         ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
% %         tremorz(1:size(ts1.data,3))=ts1.data;
% %         filt_z=filtfilt(b,a,tremorz);
% %         timeor=0:1/samplerate:(size(tremorx,2)-1)/samplerate;
% %         
% %         rep=12;
% %         for ss=1:2
% %             for ix=1:length(start{ss,1});
% %                 baseline3(1,ix)=mean(envelope(start{ss,1}(ix):ending{ss,1}(ix)));
% %                 baseline2(1,:)=envelope(start{ss,1}(ix):ending{ss,1}(ix));
% %                 for i=1:12
% %                     dum=baseline2(randi(5e4,1,rep));
% %                     hfs_p(ss,:,i)=dum;
% %                    clear dum
% %                 end
% %                 clear baseline2
% %             end
% %             HFS_median(ss,:)=median(baseline3);
% %             clear baseline3
% %         end
% %         
% %         hfs(in2,:,:)=HFS_median;
% %         hfs_tt(in2,:,:,:)=hfs_p;
% %         %     for i=1:1e6
% %         %         dum=baseline3(randi(5e4,1,rep));
% %         %         dum2=dum;
% %         %         p(i)=nanmedian(dum2);
% %         %     end
% %         %     nostim(ss,numb,:)=p;
% %         %
% %         
% %         clearvars -except  iii numb  hfs in2 hfs_tt
% %     end
% % end
% % clearvars -except hfs hfs_tt
% % cd('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data')
% % save ('HFS_median')



clear all
 cd ('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data');
% cd('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data')
load('3ax_all.mat')
load('HFS_median.mat')


for aa=1
    for m=1:2
        if (numel(find(ttall(aa,m,:)<0)))~=0 % all -except all amplifying subjects
            
            phase_peak=find(ttall(aa,m,:)==min(ttall(aa,m,:)));
             if kstest((tt1{aa,1}{m,1}(phase_peak,:))-squeeze(hfs_tt(aa,m,:,phase_peak)))==0
         [h(aa,m),p(aa,m)]=  ttest( (tt1{aa,1}{m,1}(phase_peak,:))',squeeze(hfs_tt(aa,m,phase_peak)));
             else
         [h1(aa,m),p1(aa,m)]=  signrank( (tt1{aa,1}{m,1}(phase_peak,:))',squeeze(hfs_tt(aa,m,:,phase_peak)));
             end
            clear phase_peak
        end
    end
end


%     
% 
%     if kstest((tt1{aa,1}{m,1}(phase_peak,:))'-squeeze(hfs_tt(aa,m,:,phase_peak)))==1
%         
%         for i=1:12
%             [p1(1,i),h1(1,i)]=signrank(a_s_al(:,i),a_ns_al(:,i));
%             a.s_ns=[h1(1) p1(1)];
%             
%             [p2(1,i),h2(1,i)]=signrank(a_s_al(:,i),a_s2(:,i));
%             a.s_180=[h2(1) p2(1)];
%             
%             [p3(1,i),h3(1,i)]=signrank(a_ns_al(:,i),a_ns2(:,i));
%             a.ns_180=[h3(1) p3(1)];
%             
%             
%             
%         end
%         test_a='wilcoxon';
%         
%         
%     else
%         
%         [p1,h1]=ttest(a_s_al,a_ns_al);
%         a.s_ns=[p1(1) h1(1)];
%         
%         [p2,h2]=ttest(a_s_al,a_s2);
%         a.s_180=[p2(1) h2(1)];
%         
%         [p3,h3]=ttest(a_ns_al,a_ns2);
%         a.ns_180=[p3(1) h3(1)];
%         test_a='ttest';
%     end
% 
