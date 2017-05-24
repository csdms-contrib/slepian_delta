function varargout=gulfofalaskaN(res,buf)  
% XY=GULFOFALASKAN(res,buf)
% GULFOFALASKAN(...) % Only makes a plot
%
% Returns the coordinates of glaciers in the 
% northern part of the Gulf of Alaska 
%
% INPUT:
%
% res      0 The standard, default values
%          N Splined values at N times the resolution
% buf      The region buffer you want
%
% OUTPUT:
%
% XY       Closed-curved coordinates of the continent
%
% Last modified by charig at princeton.edu, 04/04/2014

defval('res',0)
defval('buf',0)

if ~isstr(res)

% The directory where you keep the coordinates
whereitsat=fullfile(getenv('IFILES'),'GLACIERS','RGI_3_2');

% Revert to original name if unbuffered
if res==0 && buf==0
  fnpl=fullfile(whereitsat,'gulfofalaskaN.mat');
elseif buf==0;
  fnpl=fullfile(whereitsat,sprintf('%s-%i.mat','gulfofalaskaN',res));
elseif buf~=0
  fnpl=fullfile(whereitsat,sprintf('%s-%i-%g.mat','gulfofalaskaN',res,buf));
end

% If you already have a file
if exist(fnpl,'file')==2 
  load(fnpl)
  if nargout==0
    %plot(XY(:,1),XY(:,2),'k-'); axis equal; grid on
  else
    varns={XY};
    varargout=varns(1:nargout);
  end
else
  % You are about to make a file
  if res==0
    % First part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clf
    load(fullfile(getenv('IFILES'),'Glaciers','RGI_3_2','gulfofalaskaN.mat'));
    % We've already checked this file when we made it that it is correct
 
    %plot(XY(:,1),XY(:,2),'LineW',2,'Color','k');
       
  else
    XY=gulfofalaskaN(0);
    XYb=bezier(XY,res);
    XY=XYb;
  end
  
  % Do we buffer?
   if buf ~= 0
    if buf > 0
      inout='out';
    else
      inout='in';
    end
    % Make some buffered coordinates and save them for later use
    disp('Buffering the coastlines... this may take a while');
    
    % Now we can buffer this. Note that BUFFERM has gone through quite 
    % a few revisions. The cell output is no longer supported these days
    [LatB,LonB] = bufferm(XY(:,2),XY(:,1),abs(buf),inout);
    
    [latcells,loncells]=polysplit(LatB,LonB);
    
    for p = 1:length(loncells)
        theareas(p) = spharea([loncells{p} latcells{p}]);
    end
    LonB = loncells{theareas==max(theareas)};
    LatB = latcells{theareas==max(theareas)};
    
    % Periodize our way
    LonB(LonB<0) = LonB(LonB<0)+360;
    XY = [LonB LatB];
    
    % Get the unbuffered version of the neighboring region
    XY2=gulfofalaskaS(0);
    
    % Now subtract the regions
    [x,y] = polybool('subtraction',XY(:,1),XY(:,2),XY2(:,1),XY2(:,2));
    
    % A figure for test 
    figure
    plot(x,y)
    axis equal
    hold on
    plot(XY2(:,1),XY2(:,2),'r')
    
    keyboard

    % Now we look at the new piece, and we know we must fix some edges
    hdl1=figure;
    plot(x,y);
    title('This plot is used to edit the coastlines.')
    

    disp(['The functions GULFOFALASKAN has paused, and made a plot'...
    ' of the current coastlines.  These should have some artifacts that '...
    'you want to remove.  Here are the instructions to do that:'])

    disp(['DIRECTIONS:  Select the data points you want to remove with '...
        'the brush tool.  Then right click and remove them.  After you have'...
        ' finished removing the points you want, select the entire curve '...
        'with the brush tool, and type return.  The program will save the '...
        'currently brushed data in a variable, and then make another plot '...
        'for you to confirm you did it right.'])
   
    keyboard
    
    % Get the brushed data from the plot
    pause(0.1);
    hBrushLine = findall(hdl1,'tag','Brushing');
    brushedData = get(hBrushLine, {'Xdata','Ydata'});
    brushedIdx = ~isnan(brushedData{1});
    brushedXData = brushedData{1}(brushedIdx);
    brushedYData = brushedData{2}(brushedIdx);
    
    % Special case: check if the curve is closed
    if brushedXData(1)~=brushedXData(end) && brushedYData(1)~=brushedYData(end)
        brushedXData = [brushedXData brushedXData(1)];
        brushedYData = [brushedYData brushedYData(1)];
    end
    
    figure
    plot(brushedXData,brushedYData)
    title('This figure confirms the new data you selected with the brush.')
    
    disp(['The newest figure shows the data you selected with the brush '...
        'tool after you finished editing.  If this is correct, type return.'...
        '  If this is incorrect, type dbquit and run this program again to redo.'])
    keyboard
    
    XY = [brushedXData' brushedYData'];
  
   end
  
   
  % Save the file
  save(fnpl,'XY')
  
  varns={XY};
  varargout=varns(1:nargout);
  
end
  
elseif strcmp(res,'demo1')
      path(path,'~/src/m_map');
      XY1 = gulfofalaskaN(10);
      XY2 = gulfofalaskaN(10,0.5);
      figure
      m_proj('oblique mercator','longitudes',[220 220],'latitudes',[75 50],'aspect',1.0);
      m_grid;
      m_coast('color','k');
      % Original
      m_line(XY1(:,1),XY1(:,2),'color','magenta','linestyle','-');
      % Buffered
      m_line(XY2(:,1),XY2(:,2),'color','blue','linestyle','-');
      
elseif strcmp(res,'demo2')
      path(path,'~/src/m_map');
      XY1 = gulfofalaskaN(10);
      XY2 = gulfofalaskaS(10);
      figure
      m_proj('oblique mercator','longitudes',[220 220],'latitudes',[75 50],'aspect',1.0);
      m_grid;
      m_coast('color','k');
      % Original
      m_line(XY1(:,1),XY1(:,2),'color','magenta','linestyle','-');
      % Buffered
      m_line(XY2(:,1),XY2(:,2),'color','blue','linestyle','-');
  
end
  
  
end
