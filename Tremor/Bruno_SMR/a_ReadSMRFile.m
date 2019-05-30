function varargout = SmrFileRead(varargin)
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @SmrFileRead_OpeningFcn, ...
                       'gui_OutputFcn',  @SmrFileRead_OutputFcn, ...
                       'gui_LayoutFcn',  [] , ...
                    'gui_Callback',   []);
    if nargin && ischar(varargin{1})
        gui_State.gui_Callback = str2func(varargin{1});
    end

    if nargout
        [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
    else
        gui_mainfcn(gui_State, varargin{:});
    end
% End initialization code - DO NOT EDIT

% --- Executes just before SmrFileRead is made visible.
function SmrFileRead_OpeningFcn(hObject, eventdata, handles, varargin)
global defdatapn;
    handles.output = hObject;
    guidata(hObject, handles);

    defdatapn=' ';
    
% --- Outputs from this function are returned to the command line.
function varargout = SmrFileRead_OutputFcn(hObject, eventdata, handles) 
    varargout{1} = handles.output;
    
% --- Executes on button press in OpenFile.
function OpenFile_Callback(hObject, eventdata, handles)
global fid datafn defdatapn; 
       
    titText={' '};
    set(handles.graptext,'String',titText);
    
    fclose('all');
    [datafn, datapn] = uigetfile('*.smr','Select a SMR File',defdatapn);
    if isequal(datafn,0) return; end;
    defdatapn=sprintf('%s',datapn);
    FullFileName=fullfile(datapn,datafn);
    fid=fopen(FullFileName); % reading only % fid=fopen(datafn,'rb+');
    if fid > 0
      set(handles.listbox1, 'Enable', 'on');
      set(handles.listbox2, 'Enable', 'on');
      set(handles.RunBut  , 'Enable', 'on');
    else
      msgbox('Unable to open data file','File Open Error','error');
      return;
    end;
    set(handles.listbox1,'Value',1);
    set(handles.listbox2,'Value',1);
    ChanList=SONChanList(fid);
    d=size(ChanList);
    NumOfChans=d(2);
    nEvents=0; nWaves=0;
    clear ListEventsString;
    clear ListWaveString;
    for i=1:NumOfChans
        ck=ChanList(i).kind;
        switch ck
            case {1,9}      % wave data
                nWaves=nWaves+1;
                stx=sprintf('%i  %s',ChanList(i).number,ChanList(i).title);
                ListWaveString{nWaves,1}=stx;
            case {2,3,6,8}    % event,marker,
                nEvents=nEvents+1; 
                stx=sprintf('%i  %s',ChanList(i).number,ChanList(i).title);
                ListEventString{nEvents,1}=stx;
        end
    end
    
    set(handles.listbox1, 'String',ListEventString);
    set(handles.listbox2, 'String',ListWaveString);
    nString1=get(handles.listbox1, 'String');
    nString2=get(handles.listbox2, 'String');
    titCh1=nString1{get(handles.listbox1,'Value')};
    titCh2=nString2{get(handles.listbox2,'Value')};
    titText=sprintf('%s;  Ch1_%s    Ch2_%s',datafn,titCh1,titCh2);
    set(handles.graptext,'String',titText);
    
% --- Executes on button press in RunBut.

function RunBut_Callback(hObject, eventdata, handles)
global fid datafn defdatapn ; 

    set(handles.RunBut  ,'Enable', 'off');
    set(handles.OpenFile,'Enable', 'off');
    set(handles.figure1 ,'Pointer','watch');
    drawnow;
    
    nString1=get(handles.listbox1, 'String');
    nString2=get(handles.listbox2, 'String');
    titCh1=nString1{get(handles.listbox1,'Value')};
    titCh2=nString2{get(handles.listbox2,'Value')};

    [st1 st2]=strtok(titCh1,' ');
    nCh1=str2num(st1);  tCh1=st2;
    [st1 st2]=strtok(titCh2,' ');
    nCh2=str2num(st1);  tCh2=st2;

   
    %Reading Event Data
    ChanInf= SONChannelInfo(fid,nCh1);
    [data,header]=SONGetChannel(fid,nCh1);
    switch ChanInf.kind
        case {2,3} 
            MasEvent=data;
        case 6     
            MasEvent=data.timings;
    end
    
    %Reading Wave Data
    [data,header]=SONGetChannel(fid,nCh2);
    MasData=double(data);
    clear data;
    dtTime=header.sampleinterval;
    SamleRate=1.0/dtTime;
    evinfo.dtTime=dtTime;
   
    [nEvent n]=size(MasEvent);
    [n  nData]=size(MasData);
    
    tWidth=0.4;
    tOffset=0.2;
    tShiftCondition=0.000;  %10 msk
    editstr = get(handles.edit1,'string');
    timshift = str2double(editstr);
    if timshift ~=nan; tShiftCondition=timshift
%    if timshift > 0 && timshift < 1 ; tShiftCondition=timshift
    else warning('tShiftCondition was setup = 0.000'); 
    end;
      
    
    nWidth=int32(tWidth/dtTime);
    nOffset=int32(tOffset/dtTime);
    nShiftCondition=int32(tShiftCondition/dtTime);

    evinfo.tWidth=tWidth;
    evinfo.tOffset=tOffset;
    evinfo.nWidth=nWidth;
    evinfo.nOffset=nOffset;

    
    Xaxis=zeros(nWidth,1);
    s=-1*tOffset;  for i=1:nWidth  Xaxis(i)=s;  s=s+dtTime;  end

    % define posible namber of Events
    MasEnableEvents=int32(zeros(nEvent,1));
    nPosibleEvent=int32(0);
    for i = 1:nEvent;
        currtime=int32(MasEvent(i)/dtTime);
        n1=currtime-nOffset-nShiftCondition;
        n2=n1+nWidth-1+nShiftCondition;
        if n1 > 0 && n2 < nData
                nPosibleEvent=nPosibleEvent+1;
                MasEnableEvents(nPosibleEvent)=i;
        end
    end
    nEvent=nPosibleEvent;
    titText=sprintf('%s ;  Ch1_%s;  Ch2_%s  TShift=%5.3f',datafn,titCh1,titCh2,tShiftCondition);
    evinfo.grTitle=titText;
    
% Save All events  in case for later analysis, could be delete;    
    Yaxis=zeros(nEvent,1);
    EvData=zeros(nEvent,nWidth);
    for i = 1:nEvent;
        currtime=int32(MasEvent(MasEnableEvents(i))/dtTime);
        Yaxis(i)=i;
%        Yaxis(i)=currtime/SamleRate;
        n1=currtime-nOffset;
        n2=n1+nWidth-1;
        EvData(i,:) = MasData(n1:n2);
    end
    evave=mean(EvData);
%%    
    u=figure;
    rscreen=get(0,'ScreenSize');
    dxx=50;
    set(u,'Position',[1 rscreen(2)+dxx rscreen(3)/2 rscreen(4)-3*dxx] )
    set(u,'Name',titText);
    subplot(2,1,1); imagesc(Xaxis,Yaxis,EvData);
    lx=20;  caxis([-lx lx]);
    subplot(2,1,2); plot(Xaxis,evave); 
    lx=20;  ylim([-lx lx]);
%%    
    assignin('base','evinfo',evinfo);
    assignin('base','EvData',EvData);
    assignin('base','Xaxis',Xaxis);
    assignin('base','Yaxis',Yaxis);
    assignin('base','evave',evave);
    assignin('base','MasData',MasData);

    set(handles.RunBut  ,'Enable', 'on');
    set(handles.OpenFile,'Enable', 'on');
    set(handles.figure1 ,'Pointer','arrow');
    
% --- Executes on button press in CloseBut.
function CloseBut_Callback(hObject, eventdata, handles)
    close();
    
% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
global datafn; 
    h1=handles.listbox1;
    h2=handles.listbox2;
    nString1=get(h1, 'String');
    nString2=get(h2, 'String');
    titText=sprintf('%s;  Ch1:%s   Ch2:%s',datafn,nString1{get(h1,'Value')},nString2{get(h2,'Value')});
    set(handles.graptext,'String',titText);
% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
global datafn;
    h1=handles.listbox1;
    h2=handles.listbox2;
    nString1=get(h1, 'String');
    nString2=get(h2, 'String');
    titText=sprintf('%s;  Ch1:%s   Ch2:%s',datafn,nString1{get(h1,'Value')},nString2{get(h2,'Value')});
    set(handles.graptext,'String',titText);
% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global tShiftCondition;
    tShiftCondition=0.000;
    txt=sprintf('%5.3f',tShiftCondition);
    set(hObject,'String',txt);
function [c0,c1,c2,d1,d2]=a_NotchParam(sampRate,NotchFreq,epsi)
  w0=tan(pi*NotchFreq/sampRate);
  x1=w0*w0; x2=1+epsi*w0; x3=1-epsi*w0; x4=(x2*x2+x1);
  c0= (1+x1)/x4;
  c1= -2.00*(1-x1)/x4;
  c2= c0;
  d1= 2.00*(1-x1*epsi*epsi-x1)/x4;
  d2= -(x3*x3+x1)/x4;
