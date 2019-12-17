function ret=spadraw(t,x,nt,e,lnp,le,ln,it,mode,rxyz,filenm);
%SPADRAW    Draws a picture of the mechanism on the figure screen.
%
%   In simulations, only the beam and truss elements are shown
%   without numbers.
%
%   SPADRAW('filename') draws a picture using the data stored in the
%   specified SBD file. Do not specify the extension '.SBD'. When no
%   filename is specified, a file open dialog is shown.
%
%   SPADRAW() called with more paramaters is only intended for internal 
%   use e.g. by SPACAR when MODE>0 and SPASIM with plotting enabled.

%    SPACAR Program, m-file SPADRAW BETA version
%
%    Copyright 2000-2011
%    University of Twente
%    Department of Mechanical Engineering
%    Laboratory of Mechanical Automation
%
% R. Aarts 12-11-2001
% R. Aarts 17-10-2001
% R. Aarts 31- 7-2001
% H. Super 19- 7-2000
%
% For faster simulation & plotting, not every timestep will be plotted.
% Increasing normfactor generally results in less plotted timesteps and
% therefore faster simulation.
% Change plotall=0 to plotall=1 to plot mechanism at every timestep.
%
% Modified element drawing algorith for increased drawing speed
% Modified function call from SPACAR not using temporary files anymore!
% Abort button added: resulting in spacar error message "aborted by user"
% After simulation: slider.
%
% The extended calling interface depends on the number of input arguments:
% * When called with a numeric argument (i.e. not a filename this is
%   assumed to be a new time step.
% * When called with all six parameters initialization is carried out.
% * A call with 3 parameters is used for callbacks.

% To DO: Change order, improve layout, history (?), clean up

persistent SPADRAW

plotall=0;
normfactor=150;
histlenmax=1000;

retval=[];

ni = nargin;
no = nargout;

switch ni
case 0  % Open dialog to get a file
  [filename, pathname] = uigetfile('*.sbd', 'Open Spacar Binary Data file');
  if filename~=0
    try
      [pathname,filename,fileext,filever] = ...
         fileparts(fullfile(pathname,filename));
      spadraw(fullfile(pathname,[ filename filever ]));
    catch
      errordlg(lasterr,'Error Loading SBD file','modal');
%     uiwait;
    end;
  end;
case 1  % Either a "normal" user call (parameter is a filename) or
        % the next step when called from SPACAR / SPASIM.
  if isnumeric(t)
    error('Parameter must be a filename.');
  end;
  % string input: process SBD file
  filename = t;
  if exist([filename '.sbd'])~=2
    error([ 'File "' filename '.sbd" not found.']);
  end;
  % read data from SBD file
  tstp = 1;
  mstp = getfrsbf([filename '.sbd'],'tdef');
  lnp  = getfrsbf([filename '.sbd'],'lnp');
  le   = getfrsbf([filename '.sbd'],'le');
  ln   = getfrsbf([filename '.sbd'],'ln');
  it   = getfrsbf([filename '.sbd'],'it');
  mode = getfrsbf([filename '.sbd'],'mode');
  rxyz = getfrsbf([filename '.sbd'],'rxyz');
  x    = getfrsbf([filename '.sbd'],'x',tstp);
  e    = getfrsbf([filename '.sbd'],'e',tstp);
  tb   = getfrsbf([filename '.sbd'],'t',1);
  te   = getfrsbf([filename '.sbd'],'t',mstp);
  % call "normal" interface
  spadraw(tb,x,mstp,e,lnp,le,ln,it,mode,rxyz,filename);
  % modify settings
  if SPADRAW.wid>0 & ishandle(SPADRAW.wid),
    set(SPADRAW.wid,'Enable','off');
    delete(SPADRAW.wid);
  end;
  SPADRAW.wid=0;
% set(SPADRAW.gcfid,'MenuBar','figure');
  if (SPADRAW.mstp>1), 
    ss = 1/(SPADRAW.mstp-1);
    SPADRAW.sid = build_slider(SPADRAW.gcfid,tb,te,tb,ss);
  else 
    ss = 1.0; 
  end;
  SPADRAW.plota = 1;
  retval = SPADRAW.stopsp;
case 4  % Used for succesive time steps
  SPADRAW.tstp=SPADRAW.tstp+1;
  if SPADRAW.gcfid~=0 & ~ishandle(SPADRAW.gcfid),
    SPADRAW.gcfid = 0;
  end;
  if SPADRAW.gcfid~=0,
    if (max(abs(x-SPADRAW.x_old))>SPADRAW.norm | ...
        SPADRAW.plota==1 | ...
        SPADRAW.stopsp==100 | ...
        SPADRAW.tstp==SPADRAW.mstp)
      if SPADRAW.mstp>0,
        st = ['Spacar: ',SPADRAW.filename, ...
              ' [',num2str(SPADRAW.tstp),'/',num2str(SPADRAW.mstp),']'];
        set(SPADRAW.gcfid,'Name',st);
      end;
      SPADRAW.x_old=x;
      plminn = [];
      plmaxn = [];
      [plminn,plmaxn] = dominmax(SPADRAW.ndim, x(SPADRAW.xlist), ...
                                 x(SPADRAW.ylist), x(SPADRAW.zlist));
      [SPADRAW.plmin,SPADRAW.plmax] = ...
          setaxis(SPADRAW.aid, SPADRAW.ndim, plminn, plmaxn, ...
                  SPADRAW.plmin, SPADRAW.plmax, SPADRAW.tstp, SPADRAW.mstp);
      for i=1:SPADRAW.maxel
        % 2 dim
% 09-01-2008 JPM added element type 10, Planar Pinbody
% 22-01-2009 JPM added element type 15, Planar Tube
        if ( SPADRAW.it(i)==5 | SPADRAW.it(i)==7 | SPADRAW.it(i)==10 ...
           | SPADRAW.it(i)==15 )
          nodep=SPADRAW.ln(i,1);
          if (SPADRAW.it(i)==5 | SPADRAW.it(i)==10 | SPADRAW.it(i)==15),
            nodeq=SPADRAW.ln(i,3);
          else,
            nodeq=SPADRAW.ln(i,2);
          end;
          tempx=x(SPADRAW.lnp([nodep,nodeq],1))';
          tempy=x(SPADRAW.lnp([nodep,nodeq],2))';
          set(SPADRAW.lid(i),'XData',tempx,'YData',tempy);
        end;
        % 3 dim
% 18-12-2007 JPM added element type 9, Pinbody
% 28-01-2009 JPM added element type 16, Tube
% 18-03-2011 RA added element type 17, Screw
        if ( SPADRAW.it(i)==1 | SPADRAW.it(i)==2 | SPADRAW.it(i)==9 ...
           | SPADRAW.it(i)==16 | SPADRAW.it(i)==17 )
          nodep=SPADRAW.ln(i,1);
          if SPADRAW.it(i)==1,
            nodeq=SPADRAW.ln(i,3);
          elseif SPADRAW.it(i)==2,
            nodeq=SPADRAW.ln(i,2);
          elseif SPADRAW.it(i)==9,
            nodeq=SPADRAW.ln(i,3);
          else,
            nodeq=SPADRAW.ln(i,3);
          end;
          tempx=x(SPADRAW.lnp([nodep,nodeq],1))';
          tempy=x(SPADRAW.lnp([nodep,nodeq],2))';
          tempz=x(SPADRAW.lnp([nodep,nodeq],3));
          set(SPADRAW.lid(i),'XData',tempx,'YData',tempy,'ZData',tempz);
        end;
      end;
      drawnow;
    end;
    if SPADRAW.mstp==-1 & mod(SPADRAW.tstp-1,SPADRAW.histskip)==0,
      if SPADRAW.histlen==histlenmax,
        newhistlen = histlenmax/2 - 1;
        SPADRAW.thist(1:newhistlen)   = SPADRAW.thist(1:2:2*newhistlen-1);
        SPADRAW.xhist(1:newhistlen,:) = SPADRAW.xhist(1:2:2*newhistlen-1,:);
        SPADRAW.ehist(1:newhistlen,:) = SPADRAW.ehist(1:2:2*newhistlen-1,:);
        SPADRAW.histlen  = newhistlen;
        SPADRAW.histskip = SPADRAW.histskip*2;
      end;
      if mod(SPADRAW.tstp-1,SPADRAW.histskip)==0,
        SPADRAW.histlen = SPADRAW.histlen + 1;
        SPADRAW.thist(SPADRAW.histlen) = t;
        SPADRAW.xhist(SPADRAW.histlen,:) = x;
        SPADRAW.ehist(SPADRAW.histlen,:) = e;
      end;
    end;
  end;
  retval=SPADRAW.stopsp;
case 3  % Used for callbacks and special calls
  if (nt==-9)         % CloseRequestFcn
    if SPADRAW.gcfid~=0 & ishandle(SPADRAW.gcfid),
      delete(SPADRAW.gcfid);
    end;
    munlock;
    clear SPADRAW
    return;
  elseif (nt==-3)     % menu callbacks
    if (t==1)         % plotedit
      plotedit(x);
      if ~plotedit(x,'isactive')
        set(x,'Toolbar','none');
      end;
    elseif (t==2)     % printpreview
      set_slider_off(SPADRAW.sid);
      printpreview(x);
      set_slider_on(SPADRAW.sid);
    elseif (t==3)     % print
      set_slider_off(SPADRAW.sid);
      printdlg(x);
      set_slider_on(SPADRAW.sid);
    elseif (t==4)     % axes tight
      [SPADRAW.plmin,SPADRAW.plmax] = ...
         setaxis(SPADRAW.aid, SPADRAW.ndim, SPADRAW.plmin, SPADRAW.plmax, ...
                 SPADRAW.plmin, SPADRAW.plmax, -1, SPADRAW.mstp);
    end;
  elseif (nt==-2)     % termination call from SPACAR/SPASIM
    if SPADRAW.gcfid~=0 & ~ishandle(SPADRAW.gcfid),
      SPADRAW.gcfid = 0;
    end;
    if SPADRAW.gcfid~=0,
      if SPADRAW.wid>0 & ishandle(SPADRAW.wid),
        set(SPADRAW.wid,'Enable','off');
        delete(SPADRAW.wid);
      end;
      SPADRAW.wid=0;
%     set(SPADRAW.gcfid,'MenuBar','figure');
      if SPADRAW.mstp>1 & exist([SPADRAW.filename '.sbd'])==2
        SPADRAW.mstp = getfrsbf([SPADRAW.filename '.sbd'],'tdef');
        tb = getfrsbf([SPADRAW.filename '.sbd'],'t',1);
        te = getfrsbf([SPADRAW.filename '.sbd'],'t',SPADRAW.mstp);
        if (SPADRAW.mstp>1), 
          ss = 1/(SPADRAW.mstp-1);
          SPADRAW.sid = build_slider(SPADRAW.gcfid,tb,te,te,ss);
        else 
          ss = 1.0; 
        end;
        SPADRAW.plota = 1;
        x = getfrsbf([SPADRAW.filename '.sbd'],'x',SPADRAW.mstp);
        e = getfrsbf([SPADRAW.filename '.sbd'],'e',SPADRAW.mstp);
        SPADRAW.x_old = SPADRAW.x_old'; 
        SPADRAW.tstp = SPADRAW.mstp-1;
        spadraw(te,x,-1,e);
      elseif SPADRAW.mstp<0 & SPADRAW.histlen>1
        tb = SPADRAW.thist(1);
        te = SPADRAW.thist(SPADRAW.histlen);
        if (SPADRAW.histlen>1),
          ss = 1/(SPADRAW.histlen-1);
          SPADRAW.sid = build_slider(SPADRAW.gcfid,tb,te,te,ss);
        else
          ss = 1.0;
        end;
        SPADRAW.plota = 1;
        SPADRAW.mstp = -2;
      end;
    end;
  elseif (nt==-1) % Must be "ABORT" callback
    SPADRAW.stopsp=100;
  else            % Slider callback
    tmpsid = SPADRAW.sid;
    if SPADRAW.mstp>0,
      SPADRAW.tstp = round((SPADRAW.mstp-1)* ...
                           (get(tmpsid(1),'Value')-get(tmpsid(1),'Min'))/ ...
                           (get(tmpsid(1),'Max')-get(tmpsid(1),'Min')))+1;
      set(tmpsid(2), 'Position'           ,[35  20+                ...
                        (get(tmpsid(1),'Value')-get(tmpsid(1),'Min'))/   ...
                        (get(tmpsid(1),'Max')-get(tmpsid(1),'Min'))*200  ...
                                                   25 20]        , ...
                     'String'             ,num2str(get(tmpsid(1),'Value')))
      x = getfrsbf([SPADRAW.filename '.sbd'],'x',SPADRAW.tstp);
      e = getfrsbf([SPADRAW.filename '.sbd'],'e',SPADRAW.tstp);
      t = getfrsbf([SPADRAW.filename '.sbd'],'t',SPADRAW.tstp);
      SPADRAW.tstp = SPADRAW.tstp-1;
      spadraw(t,x,-1,e);
    else,
      [tmpmin,tmpind] = min(abs(SPADRAW.thist-get(tmpsid(1),'Value')));
      set(tmpsid(2), 'Position'           ,[35  20+                ...
                        (get(tmpsid(1),'Value')-get(tmpsid(1),'Min'))/   ...
                        (get(tmpsid(1),'Max')-get(tmpsid(1),'Min'))*200  ...
                                                   25 20]        , ...
                     'String'             ,num2str(SPADRAW.thist(tmpind)))
      spadraw(SPADRAW.thist(tmpind),SPADRAW.xhist(tmpind,:),-1,...
                                    SPADRAW.ehist(tmpind,:));
    end;
  end;
case 11  % Full initialization when called from SPACAR / SPASIM.
  mlock;
  SPADRAW.tstp     = 1;
  SPADRAW.mstp     = nt;
  SPADRAW.plota    = plotall;
  SPADRAW.stopsp   = 0;
  SPADRAW.lnp      = lnp;
  SPADRAW.ln       = ln;
  SPADRAW.it       = it;
  SPADRAW.filename = filenm;
  SPADRAW.histlen  = 0;
  SPADRAW.histskip = 1;
  SPADRAW.thist    = [];
  SPADRAW.xhist    = [];
  SPADRAW.ehist    = [];
  if (SPADRAW.mstp>=0 & SPADRAW.mstp<2), SPADRAW.mstp=-1; end;
  % Compute ndim
  s=size(SPADRAW.it);
  SPADRAW.maxel=s(2);
  SPADRAW.ndim=2;
% 18-12-2007 JPM Modification for added element type 9, Pinbody
% 28-01-2009 JPM Modification for added element type 16, Pinbody
% 18-03-2011 RA Modification for added element type 17, Screw
% 25-08-2017 RA Modification for added element type 18, Beamw
  for i=1:SPADRAW.maxel,
    if (SPADRAW.it(i)<4 | SPADRAW.it(i)==9 | SPADRAW.it(i)==16 ...
        | SPADRAW.it(i)==17 | SPADRAW.it(i)==18 ) 
      SPADRAW.ndim=3;
    end;
  end;
  % (Re)open figure
  if ~isfield(SPADRAW, 'gcfid'), SPADRAW.gcfid = 0; end;
  if (SPADRAW.gcfid~=0 & ishandle(SPADRAW.gcfid) &  ...
      (strncmp(get(SPADRAW.gcfid,'Name'),'Spacar: ',8) | ...
       strncmp(get(SPADRAW.gcfid,'Name'),'Spasim: ',8))),
    set(SPADRAW.gcfid,'HandleVisibility','on');
    figure(SPADRAW.gcfid);
    delete(get(SPADRAW.gcfid,'Children'));
  else,
    SPADRAW.gcfid=figure('DoubleBuffer'   , 'on', ...
                         'IntegerHandle'  , 'off',...
                         'MenuBar'        , 'none', ...
                         'NumberTitle'    , 'off', ...
                         'Tag'            , 'SpaDraw', ...
                         'Visible'        , 'off', ...
                         'CloseRequestFcn', 'spadraw(-9,-9,-9)');
  end;
  localmenu(SPADRAW.gcfid);
  % Figure title
  if SPADRAW.mstp>0,
    st = ['Spacar: ',SPADRAW.filename, ...
          ' [',num2str(SPADRAW.tstp),'/',num2str(SPADRAW.mstp),']'];
  else,
    st = ['Spasim: ',SPADRAW.filename];
    SPADRAW.histlen    = 1;
    SPADRAW.histskip   = 1;
    SPADRAW.thist(1)   = t;
    SPADRAW.xhist(1,:) = x;
    SPADRAW.ehist(1,:) = e;
  end;
  figurePos = get(SPADRAW.gcfid, 'Position');
  figurePos(3:4) = [430 300];
  set(SPADRAW.gcfid, 'Name'           , st, ...
                     'Position'       , figurePos, ...
                     'Color'          , [1 1 1], ...
                     'Visible'        , 'on');
  if SPADRAW.mstp<2,
    SPADRAW.wid=0;
  else
    SPADRAW.wid = uicontrol('Position', [10 10 50 50], ...
                            'String',   'ABORT',       ...
                            'Callback', 'spadraw(-1,-1,-1);');
  end;
% if SPADRAW.mstp==1
%   set(SPADRAW.gcfid, 'MenuBar', 'figure');
% end;
  SPADRAW.sid(1)=0;
  SPADRAW.aid = axes('DataAspectRatio', [1 1 1], ...
                     'XGrid','on', 'YGrid','on', 'ZGrid','on', ...
                     'Visible','off');
  SPADRAW.plmin = [];
  SPADRAW.plmax = [];
  SPADRAW.xlist = lnp(find(lnp(:,2)~=0 & lnp(:,4)==0),1);
  SPADRAW.ylist = lnp(find(lnp(:,2)~=0 & lnp(:,4)==0),2);
  SPADRAW.zlist = lnp(find(lnp(:,2)~=0 & lnp(:,3)~=0 & lnp(:,4)==0),3);
  [SPADRAW.plmin,SPADRAW.plmax] = ...
     dominmax(SPADRAW.ndim, x(SPADRAW.xlist), ...
              x(SPADRAW.ylist), x(SPADRAW.zlist));
  if isempty(SPADRAW.plmin), SPADRAW.plmin = zeros(1,SPADRAW.ndim); end;
  if isempty(SPADRAW.plmax), SPADRAW.plmax = zeros(1,SPADRAW.ndim); end;
  [SPADRAW.plmin,SPADRAW.plmax] = ...
     setaxis(SPADRAW.aid, SPADRAW.ndim, SPADRAW.plmin, SPADRAW.plmax, ...
             SPADRAW.plmin, SPADRAW.plmax, -1, SPADRAW.mstp);
  set(SPADRAW.aid, 'Visible','on');
  set(SPADRAW.aid, 'Position',[0.20 0.08 0.78 0.90]);

  SPADRAW.x_old = x;
  if SPADRAW.ndim==2,
    rotate3d(SPADRAW.gcfid,'off');
    for i=1:SPADRAW.maxel
% 09-01-2008 JPM added element type 10, Planar Pinbody, with colour red 
% 22-01-2009 JPM added element type 15, Planar Tube, with colour orange 
      if (it(i)==5 | it(i)==7 | it(i)==10 | it(i)==15)
        nodep = SPADRAW.ln(i,1);
        if SPADRAW.it(i)==5,
          nodeq = SPADRAW.ln(i,3);
          col='blue';
        elseif SPADRAW.it(i)==7,
          nodeq = SPADRAW.ln(i,2);
          col='green';
        elseif SPADRAW.it(i)==10,
          nodeq=SPADRAW.ln(i,3);
          col='red';
        else, 
          nodeq=SPADRAW.ln(i,3);
          col=[1.0 0.5 0.0];
        end;
        tempx = x(SPADRAW.lnp([nodep,nodeq],1))';
        tempy = x(SPADRAW.lnp([nodep,nodeq],2))';
        SPADRAW.lid(i) = line('XData',tempx,'YData',tempy, ...
                              'Color',col, 'LineWidth',2);
      end;
    end;
  else,
    set(SPADRAW.aid, 'View',[130 40]);
    RotateHndl = findobj(SPADRAW.gcfid, 'Label', '&Rotate 3D');
    set(RotateHndl, 'Enable', 'on');
% 18-12-2007 JPM added element type 9, Pinbody, with colour red 
% 28-01-2009 JPM added element type 16, Tube, with colour orange 
% 18-03-2011 RA added element type 17, Screw
% 25-08-2017 RA added element type 18, Beamw
    for i=1:SPADRAW.maxel
      if (it(i)==1 | it(i)==2 | it(i)==9 | it(i)==16 | it(i)==17 | it(i)==18)
        nodep = SPADRAW.ln(i,1);
        if SPADRAW.it(i)==1,      % Beam
          nodeq=SPADRAW.ln(i,3);
          col='blue';
        elseif SPADRAW.it(i)==2,  % Truss
          nodeq=SPADRAW.ln(i,2);
          col='green';
        elseif SPADRAW.it(i)==9,  % Pinbody
          nodeq=SPADRAW.ln(i,3);
          col='red';
        elseif SPADRAW.it(i)==16, % Tube
          nodeq=SPADRAW.ln(i,3);
          col=[1.0 0.5 0.0];
        elseif SPADRAW.it(i)==17, % Screw
          nodeq=SPADRAW.ln(i,3);
          col=[0.5 0.0 0.5];
        else %if SPADRAW.it(i)==18, % Beamw
          nodeq=SPADRAW.ln(i,4);
          col=[0.5 0.5 1.0];
        end;
        tempx = x(SPADRAW.lnp([nodep,nodeq],1));
        tempy = x(SPADRAW.lnp([nodep,nodeq],2));
        tempz = x(SPADRAW.lnp([nodep,nodeq],3));
        SPADRAW.lid(i) = line('XData',tempx,'YData',tempy,'ZData',tempz, ...
                              'Color',col, 'LineWidth',2);
      end;
    end;
  end;
  drawnow;
  SPADRAW.norm = max(SPADRAW.plmax-SPADRAW.plmin)/normfactor;
  % Hide window for "normal" plot commands. This must be called at the end!
  set(SPADRAW.gcfid, 'HandleVisibility','off');
  retval = SPADRAW.stopsp;
otherwise
  error('Wrong number of input arguments.');
end;
% return value
if no,
   ret = retval;
end;
% end spadraw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sid = build_slider(gcfid,tb,te,tcur,ss);
sid(1)=uicontrol(gcfid   , ...
                 'Position'       ,[10  10 20 250]       , ...
                 'Style'          ,'slider'              , ...
                 'String'         ,'time'                , ...
                 'Min'            ,tb                    , ...
                 'Max'            ,te                    , ...
                 'Value'          ,tcur                  , ...
                 'SliderStep'     ,[ss 10*ss]            , ...
                 'Callback'       ,'spadraw(1,1,1);'       ...
                );
sid(2)=uicontrol(gcfid   , ...
                 'Position'           ,[35  20+                ...
                    (get(sid(1),'Value')-get(sid(1),'Min'))/   ...
                    (get(sid(1),'Max')-get(sid(1),'Min'))*200  ...
                                               25 20]        , ...
                 'HorizontalAlignment','left'                , ...
                 'BackGroundColor'    ,'white'               , ...
                 'Enable'             ,'on'                  , ...
                 'String'             ,num2str(get(sid(1),'Value')) , ...
                 'Style'              ,'text'                , ...
                 'Tag'                ,'time t'                   ...
                );
sid(3)=uicontrol(gcfid   , ...
                 'Position'           ,[10 260 40 20]        , ...
                 'HorizontalAlignment','left'                , ...
                 'BackGroundColor'    ,'white'               , ...
                 'Enable'             ,'on'                  , ...
                 'String'             ,'Time'                , ...
                 'Style'              ,'text'                , ...
                 'Tag'                ,'time'                  ...
                );
sid(4)=uicontrol(gcfid   , ...
                 'Position'           ,[35   7 25 20]        , ...
                 'HorizontalAlignment','left'                , ...
                 'BackGroundColor'    ,'white'               , ...
                 'Enable'             ,'on'                  , ...
                 'String'             ,num2str(tb)           , ...
                 'Style'              ,'text'                , ...
                 'Tag'                ,'time tb'                  ...
                );
sid(5)=uicontrol(gcfid   , ...
                 'Position'           ,[35 240 25 20]        , ...
                 'HorizontalAlignment','left'                , ...
                 'BackGroundColor'    ,'white'               , ...
                 'Enable'             ,'on'                  , ...
                 'String'             ,num2str(te)           , ...
                 'Style'              ,'text'                , ...
                 'Tag'                ,'time te'                  ...
                );
% end build_slider
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [plmin,plmax]=dominmax(ndim,tempx,tempy,tempz);
if isempty(tempx),
  plmin(1) = 0;
  plmax(1) = 0;
else,
  plmin(1) = min(tempx);
  plmax(1) = max(tempx);
end;
if isempty(tempy),
  plmin(2) = 0;
  plmax(2) = 0;
else,
  plmin(2) = min(tempy);
  plmax(2) = max(tempy);
end;
if ndim==3,
  if isempty(tempz),
    plmin(3) = 0;
    plmax(3) = 0;
  else,
    plmin(3) = min(tempz);
    plmax(3) = max(tempz);
  end;
end;
%end dominmax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [plmino,plmaxo]=setaxis(aid,ndim,plminn,plmaxn,plmin,plmax,tstp,mstp);
persistent tstp_old
limstr=['Xlim';'YLim';'Zlim'];
if tstp>=0,
  if mstp>1 & tstp~=tstp_old,
    mul = max(0.1*(mstp-1)/abs(tstp-tstp_old),1);
  else,
    mul = 0;
  end;
end;
for i=1:ndim
  if plmaxn(i)-plminn(i)<3e-14,
    plminn(i) = plminn(i)-1e-14;
    plmaxn(i) = plmaxn(i)+1e-14;
  end;
  if tstp<0,
    plmino(i) = plminn(i);
    plmaxo(i) = plmaxn(i);
    set(aid,limstr(i,:),[plmino(i) plmaxo(i)]);
  else,
    plmino(i) = min(plminn(i),plmin(i));
    plmaxo(i) = max(plmaxn(i),plmax(i));
    curlim = get(aid,limstr(i,:));
    plminset = min(plmino(i),curlim(1));
    plmaxset = max(plmaxo(i),curlim(2));
    if mul>0,
      if plminset<curlim(1), plminset = mul*plmino(i)-(mul-1)*plmin(i); end;
      if plmaxset>curlim(2), plmaxset = mul*plmaxo(i)-(mul-1)*plmax(i); end;
    end;
    set(aid,limstr(i,:),[plminset plmaxset]);
  end;
end;
if tstp<0,
  tstp_old = 1;
else,
  tstp_old = tstp;
end;
%end setaxis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function localmenu(figNumber)
%  MENU Initializes our menubar.

% The Menu Labels (1:30), the CallBacks (31:50) and the Tags (51:56) for Uimenu
menustr = str2mat( ...
  '&File                                                   ', ...
  '>New                                                    ', ...
  '>Open                                                   ', ...
  '>Close                        tourbye             uimtag', ...
  '&Help                                                   ', ...
  '>MATLAB Tour Help             tour info           uimtag', ...
  '>----------------                                       ', ...
  '>MATLAB Help Window           helpwin             uimtag', ...
  '>Help Tips                    helpwin helpinfo    uimtag', ...
  '>----------------                                       ', ...
  '>Help Desk (HTML)             helpdesk            uimtag', ...
  '>Subscribe (HTML)             doc subscribe       uimtag', ...
  '&Window                                                 ', ...
  '>Main Window                  tour showmain       uimtag', ...
  '>----------------                                       ', ...
  '>Products                                         uimtag', ...
  '>>Matlab                      tour showmatlab     uimtag', ...
  '>>Toolboxes                   tour showtoolbox    uimtag', ...
  '>>Simulink                    tour showsimulink   uimtag', ...
  '>>Blocksets                   tour showblkset     uimtag', ...
  '>>Stateflow                   tour showsf               ', ...
  '>Applications                                     uimtag', ...
  '>>Analysis && Visualization   tour showdat        uimtag', ...
  '>>Mathematics                 tour showmath       uimtag', ...
  '>>Programming                 tour showprog       uimtag', ...
  '>>Dynamic System Simulation   tour showdss        uimtag', ...
  '>>Control Design              tour showcont       uimtag', ...
  '>>Signal Processing           tour showsig        uimtag', ...
  '>MATLAB in Industry           tour showindustry   uimtag', ...
  '>Resources && Contact Info    tour showcontac     uimtag'  ...
    );

menustr = str2mat( ...
  '&File                                                             ', ...
  '>Open                         spadraw                       uimtag', ...
  '>Close                        spadraw(-9,-9,-9)             uimtag', ...
  '>&Print...                    spadraw(3,gcbf,-3)            uimtag', ...
  '>Print Pre&view...            spadraw(2,gcbf,-3)            uimtag', ...
  '>Pa&ge Setup...               pagesetupdlg(gcbf)            uimtag', ...
  '>Property E&ditor...          propedit(gcbf)                uimtag', ...
  '&Tools                                                            ', ...
  '>Enable &Plot Editing         spadraw(1,gcbf,-3)            uimtag', ...
  '>&Axes tight                  spadraw(4,gcbf,-3)            uimtag', ...
  '>&Rotate 3D                   rotate3d(gcbf)                uimtag', ...
  '&Window                       winmenu(gcbo)                       ', ...
  '>blank                                                            ', ...
  '>blank                                                            ', ...
  '>blank                                                            '  ...
    );

handles = makemenu(figNumber, menustr(:,1:30), ...
                   menustr(:,31:60), menustr(:,61:66));

% Turn Enables off for some of the menus
NewHndl = findobj(handles, 'Label', 'New');
set(NewHndl, 'Enable', 'off');
RotateHndl = findobj(handles, 'Label', '&Rotate 3D');
set(RotateHndl, 'Enable', 'off');
contHndl = findobj(handles, 'Label', 'MATLAB in Industry');
set(contHndl, 'Separator', 'on');
%end localmenu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function set_slider_off(sid);
if sid(1)~=0
  set(sid,'Visible','off');
end;
%end set_slider_off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function set_slider_on(sid);
if sid(1)~=0
  set(sid,'Visible','on');
end;
%end set_slider_on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
