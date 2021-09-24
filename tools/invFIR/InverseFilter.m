function varargout = InverseFilter(varargin)


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @InverseFilter_OpeningFcn, ...
                   'gui_OutputFcn',  @InverseFilter_OutputFcn, ...
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


% --- Executes just before InverseFilter is made visible.
function InverseFilter_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

set(handles.ihL,'Visible','off');
set(handles.ihR,'Visible','off');
set(handles.ihM,'Visible','on');
set(handles.reg,'Value',1);
set(handles.ih_oct,'Value',0);
set(handles.oct,'Enable','off');

% --- Outputs from this function are returned to the command line.
function varargout = InverseFilter_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;

% --- Executes on selection change in IR.
function IR_Callback(hObject, eventdata, handles)

IR = evalin('base','who');
set(hObject,'String',IR)
list_entries = get(hObject,'String');
index_selected = get(hObject,'Value');
IR = list_entries{index_selected};
setappdata(hObject,'IR',IR);


% --- Executes during object creation, after setting all properties.
function IR_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in ih_oct.
function ih_oct_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    set(handles.oct,'Enable','on');
else
    set(handles.oct,'Enable','off');
end
 
% --- Executes on selection change in oct.
function oct_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function oct_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in win.
function win_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function win_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in ih_abs.
function ih_abs_Callback(hObject, eventdata, handles)


% --- Executes on button press in ih_cmplx.
function ih_cmplx_Callback(hObject, eventdata, handles)


function Nfft_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function Nfft_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in reg.
function reg_Callback(hObject, eventdata, handles)

if get(hObject,'Value')==1
    set(handles.f1,'Enable','on');
    set(handles.f2,'Enable','on');
    set(handles.regin,'Enable','on');
    set(handles.regout,'Enable','on');
    set(handles.minph_target,'Enable','on');
else
    set(handles.f1,'Enable','off');
    set(handles.f2,'Enable','off');
    set(handles.regin,'Enable','off');
    set(handles.regout,'Enable','off');
    set(handles.minph_target,'Enable','off');
end


function f2_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function f2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function f1_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function f1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function regout_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function regout_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function regin_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function regin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in frplot.
function frplot_Callback(hObject, eventdata, handles)
set(handles.frplot,'Value',1);
set(handles.fr2plot,'Value',0);
set(handles.irplot,'Value',0);
set(handles.ir2plot,'Value',0);
set(handles.db,'Visible','off');
plotting(handles);

% --- Executes on button press in irplot.
function irplot_Callback(hObject, eventdata, handles)
set(handles.frplot,'Value',0);
set(handles.fr2plot,'Value',0);
set(handles.irplot,'Value',1);
set(handles.ir2plot,'Value',0);
set(handles.db,'Visible','on');
plotting(handles);

% --- Executes on button press in ir2plot.
function ir2plot_Callback(hObject, eventdata, handles)
set(handles.frplot,'Value',0);
set(handles.fr2plot,'Value',0);
set(handles.irplot,'Value',0);
set(handles.ir2plot,'Value',1);
set(handles.db,'Visible','on');
plotting(handles);

% --- Executes on button press in fr2plot.
function fr2plot_Callback(hObject, eventdata, handles)
set(handles.frplot,'Value',0);
set(handles.fr2plot,'Value',1);
set(handles.irplot,'Value',0);
set(handles.ir2plot,'Value',0);
set(handles.db,'Visible','off');
plotting(handles);

% --- Executes on button press in cancel.
function cancel_Callback(hObject, eventdata, handles)
close InverseFilter

% --- Executes on button press in db.
function db_Callback(hObject, eventdata, handles)
plotting(handles);


% --- Executes on button press in close.
function close_Callback(hObject, eventdata, handles)
button = questdlg('Store filter?' );
if strcmp(button,'No')==1
    close InverseFilter
end
if strcmp(button,'Yes')==1
    pushbutton4_Callback;
end

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
ih=getappdata(handles.ok,'ih');
prompt = {'Enter name for inverste filter :'};
dlg_title = 'Name';
num_lines=1;
def = {'invfilt'};
newname = inputdlg(prompt,dlg_title,num_lines,def);
newname=cell2mat(newname);
assignin('base',newname,ih);


% --- Executes on button press in ok.
function ok_Callback(hObject, eventdata, handles)
user_entry = str2double(get(handles.L,'string'));
if isnan(user_entry) || user_entry < 0
    errordlg('Please enter a positive numeric value for filter length','Bad Input','modal')
else
if str2double(get(handles.Nfft,'string')) < str2double(get(handles.L,'string'))
     errordlg('Filter length cannot be greater than FFT length','Bad Input','modal')
else
user_entry = str2double(get(handles.Nfft,'string'));
if isnan(user_entry) || user_entry < 0
    errordlg('Please enter a positive numeric value for FFT length','Bad Input','modal')
else
    user_entry = str2double(get(handles.f1,'string'));
if isnan(user_entry) || user_entry < 0
    errordlg('Please enter a positive numeric value for f1','Bad Input','modal')
else
    user_entry = str2double(get(handles.f2,'string'));
if isnan(user_entry) || user_entry < 0
    errordlg('Please enter a positive numeric value for f2','Bad Input','modal')
else
    user_entry = str2double(get(handles.regin,'string'));
if isnan(user_entry)
    errordlg('Please enter a numeric gain value','Bad Input','modal')
else
    user_entry = str2double(get(handles.regout,'string'));
if isnan(user_entry)
    errordlg('Please enter a numeric gain value','Bad Input','modal')
else
    
    IR=getappdata(handles.IR,'IR');
    Nfft=str2double(get(handles.Nfft,'string'));
    L=str2double(get(handles.L,'string'));
    
    if get(handles.win,'Value')==1
        win=1;
    else
        win=0;
    end

    if get(handles.reg,'Value')==0
        f1=0;
        f2=22050;
        reg_in=inf;
        reg_out=inf;
    else
        f1=str2double(get(handles.f1,'string'));
        f2=str2double(get(handles.f2,'string'));
        reg_in=str2double(get(handles.regin,'string'));
        reg_out=str2double(get(handles.regout,'string'));
    end
    if get(handles.ih_oct,'Value')==1  
        Noct=str2double(get(handles.oct,'String'));
    else
        Noct=0;
    end
    if get(handles.ih_cmplx,'Value')==1
        type='complex';
    elseif get(handles.ih_abs,'Value')==1
         type='linphase';
    elseif get(handles.ih_minph,'Value')==1
         type='minphase';
    end
    
    h=evalin('base',IR);
    [ih]=invFIR(type,h,Nfft,Noct,L,[f1 f2],[reg_in reg_out],win);
    convh=ifft(fft(h,(length(h)+length(ih)-1)).*fft(ih,(length(h)+length(ih)-1)));
    setappdata(hObject,'ih',ih);
    setappdata(hObject,'convh',convh./max(max(abs(convh))));
    setappdata(hObject,'h',h);
    setappdata(handles.ok,'Nfft',Nfft);
    plotting(handles);
end
end
end
end
end
end
end


function L_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function L_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in ih_minph.
function ih_minph_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    set(handles.minph_target,'Enable','off');
    set(handles.reg,'Value',1);
    set(handles.reg,'Enable','off');
end

function plotting(handles)
cla;
Nfft=getappdata(handles.ok,'Nfft');
ih=getappdata(handles.ok,'ih');
h=getappdata(handles.ok,'h');
convh=getappdata(handles.ok,'convh');
iH=20*log10(abs(fft(ih)));
convH=20*log10(abs(fft(convh,Nfft)));
H=20*log10(abs(fft(h,Nfft)));
freq=0:44100/(Nfft-1):22050;
ih=ih./max(max(abs(ih)));

if get(handles.db,'Value')==1
    ih=20*log10(abs(ih));
    convh=20*log10(abs(convh));
end
if size(ih,2)==1
    set(handles.ihL,'Visible','off');
    set(handles.ihR,'Visible','off');
    set(handles.ihM,'Visible','on');
else
    set(handles.ihL,'Visible','on');
    set(handles.ihR,'Visible','on');
    set(handles.ihM,'Visible','off');
end

% case plot filter impulse response
if get(handles.irplot,'Value')==1
    if size(ih,2)==2
        axes(handles.ihL);
        cla(handles.ihL,'reset');
        plot(ih(:,1),'Color',[0 0 0]);
        set(handles.ihL,'XGrid','on','FontWeight','bold');
        axis([0 length(ih) -1 1]);
        if get(handles.db,'Value')==1
            axis([0 length(ih) -60 0]);
        end
        box('on');
        hold('all');
        ylabel('CH1 amplitude','FontWeight','bold');
        
        axes(handles.ihR);
        cla(handles.ihR,'reset');
        plot(ih(:,2),'Color',[0 0 0]);
        set(handles.ihR,'XGrid','on','FontWeight','bold');
        axis([0 length(ih) -1 1]);
        if get(handles.db,'Value')==1
            axis([0 length(ih) -60 0]);
        end
        box('on');
        hold('all');
        xlabel('time / samples','FontWeight','bold');
        ylabel('CH2 amplitude','FontWeight','bold');
    else
        axes(handles.ihM);
        cla(handles.ihM,'reset');
        plot(ih,'Color',[0 0 0]);
        set(handles.ihM,'XGrid','on','FontWeight','bold');
        axis([0 length(ih) -1 1]);
        if get(handles.db,'Value')==1
            axis([0 length(ih) -60 0]);
        end
        box('on');
        hold('all');
        xlabel('time / samples','FontWeight','bold');
        ylabel('amplitude','FontWeight','bold');
    end
end
        

% case plot filter frequency response
if get(handles.frplot,'Value')==1
     if size(ih,2)==2
        axes(handles.ihL);
        cla(handles.ihL,'reset');
        set(handles.ihL,'XGrid','on','FontWeight','bold');
        set(get(gca,'YLabel'),'String','CH1 amplitude / dB','FontSize',10,'FontWeight','bold')
        set(handles.ihL,...
        'FontWeight','bold',...
        'XGrid','on',...
        'XMinorGrid','on',...
        'XTick',[32,63,125,250,500,1000,2000,4000,8000,16000],...
        'XTickLabel',{'32','63','125','250','500','1k','2k','4k','8k','16k'},...
        'XScale','log','YGrid','on')
        XLim([32 16000]);
        box('on');
        hold('all');
        plot(freq,-H(1:end/2,1),'Color',[1 0 0],'DisplayName','true inverse');
        plot(freq,iH(1:end/2,1),'Color',[0 0 1],'DisplayName','designed inverse');
        legend1 = legend('show');
        set(legend1,'Location','SouthWest','FontSize',8);
        
        axes(handles.ihR);
        cla(handles.ihR,'reset');
        set(handles.ihR,'XGrid','on','FontWeight','bold');
        set(get(gca,'YLabel'),'String','CH2 amplitude / dB','FontSize',10,'FontWeight','bold')
        set(get(gca,'XLabel'),'String','frequency / Hz','FontSize',10,'FontWeight','bold')
        set(handles.ihR,...
        'FontWeight','bold',...
        'XGrid','on',...
        'XMinorGrid','on',...
        'XTick',[32,63,125,250,500,1000,2000,4000,8000,16000],...
        'XTickLabel',{'32','63','125','250','500','1k','2k','4k','8k','16k'},...
        'XScale','log','YGrid','on')
        XLim([32 16000]);
        box('on');
        hold('all');
        plot(freq,-H(1:end/2,2),'Color',[1 0 0],'DisplayName','true inverse');
        plot(freq,iH(1:end/2,2),'Color',[0 0 1],'DisplayName','designed inverse');
        legend1 = legend('show');
        set(legend1,'Location','SouthWest','FontSize',8);
    else
        axes(handles.ihM);
        cla(handles.ihM,'reset');
        set(handles.ihM,'XGrid','on','FontWeight','bold');
        set(get(gca,'YLabel'),'String','amplitude / dB','FontSize',10,'FontWeight','bold')
        set(get(gca,'XLabel'),'String','frequency / Hz','FontSize',10,'FontWeight','bold')
        set(handles.ihM,...
        'FontWeight','bold',...
        'XGrid','on',...
        'XMinorGrid','on',...
        'XTick',[32,63,125,250,500,1000,2000,4000,8000,16000],...
        'XTickLabel',{'32','63','125','250','500','1k','2k','4k','8k','16k'},...
        'XScale','log','YGrid','on')
        XLim([32 16000]);
        box('on');
        hold('all');
        plot(freq,-H(1:end/2),'Color',[1 0 0],'DisplayName','true inverse');
        plot(freq,iH(1:end/2),'Color',[0 0 1],'DisplayName','designed inverse');
        legend1 = legend('show');
        set(legend1,'Location','SouthWest','FontSize',8);
    end
end


% case plot filtered impulse response
if get(handles.ir2plot,'Value')==1
    if size(convh,2)==2
        axes(handles.ihL);
        cla(handles.ihL,'reset');
        plot(convh(:,1),'Color',[0 0 0]);
        set(handles.ihL,'XGrid','on','FontWeight','bold');
        axis([0 length(convh) -1 1]);
        if get(handles.db,'Value')==1
            axis([0 length(convh) -60 0]);
        end
        box('on');
        hold('all');
        ylabel('CH1 amplitude','FontWeight','bold');
        
        axes(handles.ihR);
        cla(handles.ihR,'reset');
        plot(convh(:,2),'Color',[0 0 0]);
        set(gca,'XGrid','on','FontWeight','bold');
        axis([0 length(convh) -1 1]);
        if get(handles.db,'Value')==1
            axis([0 length(convh) -60 0]);
        end
        box('on');
        hold('all');
        xlabel('time / samples','FontWeight','bold');
        ylabel('CH2 amplitude','FontWeight','bold');
    else
        axes(handles.ihM);
        cla(handles.ihM,'reset');
        plot(convh,'Color',[0 0 0]);
        set(gca,'XGrid','on','FontWeight','bold');
        axis([0 length(convh) -1 1]);
        if get(handles.db,'Value')==1
            axis([0 length(convh) -60 0]);
        end
        box('on');
        hold('all');
        xlabel('time / samples','FontWeight','bold');
        ylabel('amplitude','FontWeight','bold');
    end
end
%     
% case plot filtered frequency response
if get(handles.fr2plot,'Value')==1
     if size(ih,2)==2
        axes(handles.ihL);
        cla(handles.ihL,'reset');
        set(handles.ihL,'XGrid','on','FontWeight','bold');
        set(get(gca,'YLabel'),'String','CH1 amplitude / dB','FontSize',10,'FontWeight','bold')
        set(handles.ihL,...
        'FontWeight','bold',...
        'XGrid','on',...
        'XMinorGrid','on',...
        'XTick',[32,63,125,250,500,1000,2000,4000,8000,16000],...
        'XTickLabel',{'32','63','125','250','500','1k','2k','4k','8k','16k'},...
        'XScale','log','YGrid','on')
        XLim([32 16000]);
        box('on');
        hold('all');
        plot(freq,H(1:end/2,1),'Color',[1 0 0],'DisplayName','original');
        plot(freq,convH(1:end/2,1),'Color',[0 0 1],'DisplayName','filtered');
        legend1 = legend('show');
        set(legend1,'Location','SouthWest','FontSize',8);
        
        axes(handles.ihR);
        cla(handles.ihR,'reset');
        set(handles.ihR,'XGrid','on','FontWeight','bold');
        set(get(gca,'YLabel'),'String','CH2 amplitude / dB','FontSize',10,'FontWeight','bold')
        set(get(gca,'XLabel'),'String','frequency / Hz','FontSize',10,'FontWeight','bold')
        set(handles.ihR,...
        'FontWeight','bold',...
        'XGrid','on',...
        'XMinorGrid','on',...
        'XTick',[32,63,125,250,500,1000,2000,4000,8000,16000],...
        'XTickLabel',{'32','63','125','250','500','1k','2k','4k','8k','16k'},...
        'XScale','log','YGrid','on')
        XLim([32 16000]);
        box('on');
        hold('all');
        plot(freq,H(1:end/2,2),'Color',[1 0 0],'DisplayName','original');
        plot(freq,convH(1:end/2,2),'Color',[0 0 1],'DisplayName','filtered');
        legend1 = legend('show');
        set(legend1,'Location','SouthWest','FontSize',8);
    else
        axes(handles.ihM);
        cla(handles.ihM,'reset');
        set(handles.ihM,'XGrid','on','FontWeight','bold');
        set(get(gca,'YLabel'),'String','amplitude / dB','FontSize',10,'FontWeight','bold')
        set(get(gca,'XLabel'),'String','frequency / Hz','FontSize',10,'FontWeight','bold')
        set(handles.ihM,...
        'FontWeight','bold',...
        'XGrid','on',...
        'XMinorGrid','on',...
        'XTick',[32,63,125,250,500,1000,2000,4000,8000,16000],...
        'XTickLabel',{'32','63','125','250','500','1k','2k','4k','8k','16k'},...
        'XScale','log','YGrid','on')
        XLim([32 16000]);
        hold('all');
        plot(freq,H(1:end/2),'Color',[1 0 0],'DisplayName','original');
        plot(freq,convH(1:end/2),'Color',[0 0 1],'DisplayName','filtered');
        legend1 = legend('show');
        set(legend1,'Location','SouthWest','FontSize',8);
     end
end
