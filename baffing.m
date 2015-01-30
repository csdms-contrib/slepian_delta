function varargout=baffin(res,buf)  
% XY=BAFFIN(res,buf)
% 
% Finds the coordinates of Baffin Island, taking into account
% the fact that it is next to Ellesmere Island, which 
% causes problems when you buffer.
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
% NOTES: We start with a region outline created by encircling glaciers in
% the Randolph Glacier Inventory 3.2.  This outline is buffered as
% requested.  Since buffered regions infringe on neighboring regions, we
% subtract the outlines for Ellesmere Island.  This can cause
% some artifacts near some of the boundary intersections, so this function
% employs the same brushing technique as the basin functions for
% Antarctica.  We make the region, plot it, remove erroneous points using
% the brushing tool, and save this new outline.
%
% Last modified by charig at princeton.edu, 09/04/2014

defval('res',10)
defval('buf',0)

if ~isstr(res) % Not a demo

  % The directory where you keep the coordinates
  whereitsat=fullfile(getenv('IFILES'),'GLACIERS','RGI_3_2');

  % Revert to original name if unbuffered
  if res==0 && buf==0
    fnpl=fullfile(whereitsat,'Baffin.mat');
  elseif buf==0;
    fnpl=fullfile(whereitsat,sprintf('%s-%i.mat','Baffin',res));
  elseif buf~=0
    fnpl=fullfile(whereitsat,sprintf('%s-%i-%g.mat','Baffin',res,buf));
  end

  % If you already have a file
  if exist(fnpl,'file')==2 
    load(fnpl)
    if nargout==0
      plot(XY(:,1),XY(:,2),'k-'); axis equal; grid on
    else
      varns={XY};
      varargout=varns(1:nargout);
    end
  else
    % You are about to make a file
  
    % Do we buffer? Not here, so do it regular
    if buf==0
        if res==0
          % First part 
          load(fullfile(getenv('IFILES'),'Glaciers','RGI_3_2','Baffin.mat'));
          % We've already checked this file when we made it that it is correct   
        else
          XY=baffin(0);
          XYb=bezier(XY,res);
          XY=XYb;
        end
    end
  
    if buf ~=0
      % Check if we have this buffer already but only at base resolution, 
      % and then change the res on that file
      fnpl2=fullfile(whereitsat,sprintf('%s-%i-%g.mat','Baffin',0,buf));
      if exist(fnpl2,'file')==2
        load(fnpl2)
        XYb=bezier(XY,res);
        XY=XYb;
      else
        % We make a new buffer
        disp('Buffering the coastlines... this may take a while');
        XY=baffin(10);
        if buf > 0
           inout='out';
        else
           inout='in';
        end
      
        %%%
        % First we get the thing next door (Ellesmere)
        %%%
        XY1 = baffin(10);
        [LatB,LonB] = bufferm(XY1(:,2),XY1(:,1),0.2,'out');
        % Note that, if due to BEZIER there might be a pinched-off loop in
        % the XY you will get an extra NaN and will need to watch it
        % If this shouldn't work, keep it all unchanged in other words
        try
          % You'll need a line for every possible version behavior
          % Note how POLY2CW has disappeared from BUFFERM
          if sign(buf)<0 || ~isempty(strfind(version,'2010a')) 
	        % Take the last bit of non-NaNs; there might have been pinched
            % off sections in the original
	        LonB=LonB(indeks(find(isnan(LonB)),'end')+1:end);
	        LatB=LatB(indeks(find(isnan(LatB)),'end')+1:end);
          elseif ~isempty(strfind(version,'2011a')) || ~isempty(strfind(version,'2012a'))
	        LonB=LonB(1:find(isnan(LonB))-1);
	        LatB=LatB(1:find(isnan(LatB))-1);
          end
        catch
          disp('BUFFERM failed to buffer as expected')
        end
        % Now subtract a smaller Ellesmere from Greenland
        % Note: Say you ask for buf=1.0.  Then we are making Ellesmere1.0
        % and we subtract from it Baffin0.2.  This is considered the
        % thing next door.  If we just subtracted a Ellesmere1.0, we
        % would lose large sections of Baffin.
        LonB(LonB<0) = LonB(LonB<0)+360;
        XY1 = [LonB LatB];
        XY2 = ellesmere(10,buf);
        [x1,y1] = polybool('subtraction',XY2(:,1),XY2(:,2),XY1(:,1),XY1(:,2));
        % Do you need to see it to understand? Use the code below
        %path(path,'~/src/m_map');
        %figure
        %m_proj('oblique mercator','longitudes',[318 318],'latitudes',[90 50],'aspect',1.0);
        %m_grid;         %m_coast('color','k');
        %m_line(x1(:,1),y1(:,2),'color','magenta','linestyle','-');
      
        %%%
        % Now we get what we want (Baffin)
        %%%
        [LatB,LonB] = bufferm(XY(:,2),XY(:,1),abs(buf),inout);
        % Check for pinched loops
        try
          % You'll need a line for every possible version behavior
          % Note how POLY2CW has disappeared from BUFFERM
          if sign(buf)<0 || ~isempty(strfind(version,'2010a')) 
            % Take the last bit of non-NaNs; there might have been pinched
            % off sections in the original
            LonB=LonB(indeks(find(isnan(LonB)),'end')+1:end);
            LatB=LatB(indeks(find(isnan(LatB)),'end')+1:end);
          elseif ~isempty(strfind(version,'2011a')) || ~isempty(strfind(version,'2012a'))
            LonB=LonB(1:find(isnan(LonB))-1);
            LatB=LatB(1:find(isnan(LatB))-1);
          end
        catch
          disp('BUFFERM failed to buffer as expected')
        end
        % Periodize our way
        LonB(LonB<0) = LonB(LonB<0)+360;
        % Now subtract the thing next door so we don't overlap Ellesmere
        [x,y] = polybool('subtraction',LonB,LatB,x1,y1);  
    
        %%%
        % Assembly complete!
        %%%
 
        % Now we look at the new piece, and we know we must fix some edges
        % where the things intersected.
        hdl1=figure;
        plot(x,y);
        title('This plot is used to edit the coastlines.')
    
        fprintf(['The functions BAFFIN has paused, and made a plot \n'...
        ' of the current coastlines.  These should have some artifacts that \n'...
        'you want to remove.  Here are the instructions to do that:\n \n'])

        fprintf(['DIRECTIONS:  Select the data points you want to remove with \n'...
        'the brush tool.  Then right click and remove them.  After you have\n'...
        ' finished removing the points you want, select the entire curve \n'...
        'with the brush tool, and type return.  The program will save the \n'...
        'currently brushed data in a variable, and then make another plot \n'...
        'for you to confirm you did it right.\n'])
        keyboard
    
        % Get the brushed data from the plot
        pause(0.1);
        hBrushLine = findall(hdl1,'tag','Brushing');
        brushedData = get(hBrushLine, {'Xdata','Ydata'});
        brushedIdx = ~isnan(brushedData{1});
        brushedXData = brushedData{1}(brushedIdx);
        brushedYData = brushedData{2}(brushedIdx);
    
        figure
        plot(brushedXData,brushedYData)
        title('This figure confirms the new data you selected with the brush.')
    
        fprintf(['The newest figure shows the data you selected with the brush \n'...
        'tool after you finished editing.  If this is correct, type return.\n'...
        '  If this is incorrect, type dbquit and run this program again to redo.\n'])
        keyboard
    
        XY = [brushedXData' brushedYData'];
        
        % If we get too big, Greenland intrudes, so clip this.
        if buf > 1.5
            XY3 = greenland(10,buf,1);
            [x,y] = polybool('subtraction',XY(:,1),XY(:,2),XY3(:,1),XY3(:,2));   
            XY = [x y];
        end
  
      end % end if fnpl2 exist
    
    end % end if buf>0
   
    % Save the file
    save(fnpl,'XY')
  
    varns={XY};
    varargout=varns(1:nargout);
  
  end % end if fnpl exist
  
elseif strcmp(res,'demomake')
      % This demo illustrates the proper order that is needed to make these
      % coordinate files, since there are dependancies.
      % The coordinates need to be created in order of increasing buffer,
      % and ellesmere before baffin
      XY = ellesmere(10,0.2);
      XY = baffin(10,0.2);
      XY = ellesmere(10,0.5);
      XY = baffin(10,0.5);
      XY = ellesmere(10,1);
      XY = baffin(10,1);
      XY = ellesmere(10,1.5);
      XY = baffin(10,1.5);
      XY = ellesmere(10,2.0);
      XY = baffin(10,2.0);
  
elseif strcmp(res,'demo1')
      path(path,'~/src/m_map');
      XY1 = baffin(10);
      XY3 = ellesmere(10,2);
      XY2 = baffin(10,2);
      XY4 = greenland(10,2,1);
      figure
      m_proj('oblique mercator','longitudes',[318 318],'latitudes',[90 50],'aspect',1.0);
      m_grid;
      m_coast('color','k');
      m_line(XY3(:,1),XY3(:,2),'color','green','linestyle','-');
      XY3 = ellesmere(10,1);      
      m_line(XY3(:,1),XY3(:,2),'color','green','linestyle','-');
      % Original
      m_line(XY1(:,1),XY1(:,2),'color','magenta','linestyle','-');
      % Buffered
      m_line(XY2(:,1),XY2(:,2),'color','blue','linestyle','-');
      m_line(XY4(:,1),XY4(:,2),'color','red','linestyle','-');

      
      
elseif strcmp(res,'demo2')
      path(path,'~/src/m_map');
      XY1 = ellesmere(10,0.2);
      XY2 = ellesmere(10,0.5);
      XY3 = greenland(10,0.5);
      
      [x1,y1] = polybool('subtraction',XY3(:,1),XY3(:,2),XY1(:,1),XY1(:,2));
      [x2,y2] = polybool('subtraction',XY2(:,1),XY2(:,2),x1,y1);
      
      figure
      m_proj('oblique mercator','longitudes',[318 318],'latitudes',[90 50],'aspect',1.0);
      m_grid;
      m_coast('color','k');
      % Original
      m_line(x2,y2,'color','magenta','linestyle','-');
      % Buffered
      %m_line(XY2(:,1),XY2(:,2),'color','blue','linestyle','-');
      %m_line(XY3(:,1),XY3(:,2),'color','green','linestyle','-');
  
end
  
  
end
