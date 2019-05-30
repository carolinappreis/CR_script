    fclose('all');
    clear all;
    disp('--------------- START -');    
    defdatapn=' ';
    [datafn, datapn] = uigetfile('*.smr','Select a SMR File',defdatapn);
    if isequal(datafn,0) return; end;
    defdatapn=sprintf('%s',datapn);
    FullFileName=fullfile(datapn,datafn);
    fid=fopen(FullFileName); % reading only % fid=fopen(datafn,'rb+');
    if fid <= 0
      msgbox('Unable to open data file','File Open Error','error');
      return;
    end;
    ChanList=SONChanList(fid);
    d=size(ChanList);
    NumOfChans=d(2);
    nEvents=0; nWaves=0;
    clear EventNumbers EventTitles;
    clear WaveNumbers  WaveTitles;
    for i=1:NumOfChans
        ck=ChanList(i).kind;
        switch ck
            case {1,9}      % wave data
                nWaves=nWaves+1;
                WaveNumbers(nWaves)=ChanList(i).number;
                stx=sprintf('%s',ChanList(i).title);
                WaveTitles{nWaves,1}=stx;
            case {2,3,6}    % event,marker,
                nEvents=nEvents+1; 
                EventNumbers(nEvents)=ChanList(i).number;
                stx=sprintf('%s',ChanList(i).title);
                EventTitles{nEvents,1}=stx;
        end
    end
    
%% Reading Wave Data if there are any 
    WaveData=[];
    if nWaves > 0 
        for ix=1:nWaves
            ch=WaveNumbers(ix);
            [data,header]=SONGetChannel(fid,ch);
            l(ix)=length(data);
        end
        for ix=1:nWaves
            ch=WaveNumbers(ix);
            [data,header]=SONGetChannel(fid,ch);
            dtTime=header.sampleinterval;
            if ix == 1; SamleRate=1.0/dtTime; end;
            [nep  nData]=size(data);
            stx=sprintf('WaveChan %d %s dim=%d  SamleRate=%f',WaveNumbers(ix),WaveTitles{ix},nep,SamleRate);
            disp(stx);
            if nep ==1; WaveData(ix,:)=double(data(1:min(l)));
            else
                ln=1;
                for n=1:nep
                    nps=header.npoints(n);
                    WaveData(ix,ln:(ln+nps-1))=double(data(n,1:nps));
                    ln=ln+nps;
                end
                nData=ln;
            end  
            if header.kind == 1;  % Channel is Waveform
               WaveData(ix,:)=WaveData(ix,:)*header.scale/6553.6+header.offset;
            end;
        end;    
    end;    
%% Reading Event Data if there are any 
    clear EventData;
    if nEvents > 0 
        for ix=1:nEvents
            nEvCh=EventNumbers(ix);
            ChanInf= SONChannelInfo(fid,nEvCh);
            [data,header]=SONGetChannel(fid,nEvCh);
            switch ChanInf.kind
                case {2,3} 
                    MasEvent=data;
                case {5,6}
                    MasEvent=data.timings;
                    MasMarks=data.markers;
            end  
            EventData(ix)={MasEvent};
            nevs=length(MasEvent);
            stx=sprintf('EventChan %d  %s  nEvs=%d',nEvCh,EventTitles{ix},nevs);
            disp(stx);
        end
    end
    
%%    
    fclose('all');  
%%    
    clear SmrData;
    SmrData.FileName=FullFileName;
    SmrData.SR=SamleRate;
    SmrData.WvTits=WaveTitles;
    SmrData.WvData=WaveData;
    if nEvents > 0
        SmrData.EvTits=EventTitles;
        SmrData.EvData=EventData;
    end;
    
    SmrMatFile=strtok(FullFileName,'.');
    SmrMatFile=sprintf('%s.mat',SmrMatFile);
    save(SmrMatFile,'SmrData', '-v7.3');
    disp('--------------- DONE -');    
%%
%     [nCh,ndata]=size(SmrData.WvData);
%     Xaxis=(1:nData)/SamleRate;
%     figure;
%     for ich=1:nCh
%         subplot(nCh,1,ich);
%         plot(Xaxis,SmrData.WvData(ich,:)); xlim([Xaxis(1) Xaxis(end)]);  grid on;
%         title(WaveTitles{ich,1});
%     end;    
%     title(FullFileName)
