function WaveData= Read_SMR_File(FullFileName)
%     fclose('all');
% %     clear all;
%     defdatapn=' ';
%     [datafn, datapn] = uigetfile('*.smr','Select a SMR File',defdatapn);
%     if isequal(datafn,0) return; end;
%     defdatapn=sprintf('%s',datapn);
%     FullFileName=fullfile(datapn,datafn);
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
    %Reading Event Data
    EventData=[];
    if nEvents > 0
        for ix=1:nEvents
            ch=EventNumbers(ix);
            ChanInf= SONChannelInfo(fid,ch);
            [data,header]=SONGetChannel(fid,ch);
            switch ChanInf.kind
                case {2,3};  MasEvent=data;
                case 6    ;  MasEvent=data.timings;
            end
            EventData(ix,:)=data;
        end;   
    end;    
    
    %Reading Wave Data
    if nWaves > 0 
        for ix=1:nWaves
            ch=WaveNumbers(ix);
            [data,header]=SONGetChannel(fid,ch);
            dtTime=header.sampleinterval;
            SamleRate=1.0/dtTime;
            [nep  nData]=size(data);
            disp('__________________________________________________________________');
            stx=sprintf('Chan= %d ,  Title=  %s ,   Dim= %d ,',WaveNumbers(ix),WaveTitles{ix},nep);
            disp(stx);
            disp(header.mode);
            disp('Part    StartTime  StopTime');
            for n=1:nep
                stx=sprintf('%d  %f  %f',n,header.start(n),header.stop(n));
                disp(stx);
            end
            
            if ix==1; 
                if nep==1 
                    nDataMax=nData;
                else
                    tDataMax=max(header.stop);
                    nDataMax=tDataMax*SamleRate
                end
                WaveData=double(zeros(nWaves,nDataMax)); 
            end;
            
            if nep==1; WaveData(ix,1:length(data))=double(data);
            else
                ln=1;
                for n=1:nep
                    nps=header.npoints(n);
                    WaveData(ix,ln:(ln+nps-1))=double(data(n,1:nps));
                    ln=ln+nps;
                end
                nData=ln;
            end   
        end;    
    end;    
   
    
    disp('__________________________________________________________________');
    disp('Filal info');
    [nEvent n]=size(EventData)
    dtTime
    [n  nData]=size(WaveData);
    datasec=dtTime*nData
    nDataMax
    nData
    fclose('all');
    
    return