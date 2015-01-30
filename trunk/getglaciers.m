function varargout=getglaciers(region)
% XY=GETGLACIERS(region)
%
% Return the XY coordinates for glaciers from the Randolph Glacier
% Inventory, either for a specific region or for the entire globe.
%
% INPUT:
% 
% region    The single-region string for which you want glacier coordinates:
%           'Alaska', 'AntarcticSubantarctic', 'ArcticCanadaNorth',
%           'ArcticCanadaSouth', 'CaucasusMiddleEast', 'CentralAsia',
%           'CentralEurope', 'GreenlandPeriphery', 'Iceland',
%           'LowLatitudes', 'NewZealand', 'NorthAsia', 'RussianArctic',
%           'Scandinavia', 'SouthAsiaEast', 'SouthAsiaWest',
%           'SouthernAndes', 'Svalbard', 'WesternCanadaUS', OR:
%           'globe' in which case you get all of the above [default]
%
% OUTPUT:
%
% XY    The coordinates for the glaciers in the region you requested,
%       between 0 and 360 degrees longitude.
%     
% NOTE: 
%
% The Randolph Glacier Inventory is found at http://www.glims.org/RGI/.  
% The data are made available under certain use constraints, and users
% should consult the website for valid purposes.
% 
% Open a MATLABPOOL and this will run in parallel, automatically
%
% SEE ALSO: 
%
% SHAPEREAD from the mapping toolbox is not sufficient; instead,
% M_SHAPEREAD from the freeware M_MAP needs to be used
% (http://www2.ocgy.ubc.ca/~rich/map.html) 
%
% Last modified by charig-at-princeton.edu, 08/19/2014
% Last modified by fjsimons-at-princeton.edu, 09/23/2014

% Determine parameters and set defaults
defval('xver',0);
defval('region','globe');
defval('XYtot',[]);

% Where and under what name you find the old and save the new files
gdir=fullfile(getenv('IFILES'),'GLACIERS','RGI_3_2');
ddir1=fullfile(gdir,'SHPFILES',filesep);
ddir2=fullfile(gdir,'MATFILES',filesep);
fnpl=sprintf('%s%s.mat',ddir2,region);

% If this file already exists, load it.  Otherwise, or if we force it, make
% a new one (e.g. you added extra months to the database).
if exist(fnpl,'file')==2
  load(fnpl)
  disp(sprintf('%s loaded by GETGLACIERS',fnpl))
else
  % Make the cell of regions if you want all of them
  if strcmp(region,'globe')
    region = {'Alaska' 'ArcticCanadaNorth' 'WesternCanadaUS' 'ArcticCanadaSouth'...
	      'GreenlandPeriphery' 'Iceland' 'Svalbard' 'Scandinavia' 'RussianArctic'...
	      'NorthAsia' 'CentralEurope' 'CaucasusMiddleEast' 'CentralAsia'...
	      'SouthAsiaWest' 'SouthAsiaEast' 'LowLatitudes' 'SouthernAndes'...
	      'NewZealand' 'AntarcticSubantarctic'};
  end
  % If there is only one of them, in a specific string
  if isstr(region)
    disp(sprintf('\nGETGLACIERS making %s\n',fnpl))

    % Get the file basename (make sure it's unique so it ends in .shp)
    wholefile=ls2cell([ddir1 '*' region '.shp'],1);
    % Shaperead it! No need to use FILEPARTS, is there? Only take first
    % entry in case there are more, which there shouldn't be
    length(wholefile)
    S=m_shaperead(pref(wholefile{1}));
    % Make sure each line ends with a NaN
    coordsnan=cellfun(@(x) [x(:,1:2); NaN NaN],S.ncst,'UniformOutput',false); 
    % Now make it one vector
    X=indeks(cell2mat(coordsnan),':,1');
    Y=indeks(cell2mat(coordsnan),':,2');
    
    % See if it still works when we reduce it 
    [Y,X,cerr,tol]=reducem(Y,X);
    XY=[X Y];
    
    % Make sure the coordinates make sense
    XY(:,1)=XY(:,1)-360*[XY(:,1)>360];
    XY(:,1)=XY(:,1)+360*[XY(:,1)<0];
    % Run penlift with a very large median filter
    [xdata,ydata]=penlift(XY(:,1),XY(:,2),200);
    XY=[xdata ydata];
    
    % Save this in the file
    save(fnpl,'XY')
  elseif iscell(region)
    % Then do all of them, recursively, and in parallel!
    parfor j=1:length(region)
      XY{j}=getglaciers(region{j});
    end
    % Then stitch them all together
    % Probably best to initialize them
    lens=cellfun(@(x) length(x),XY);
    for k=1:length(region)
      XYtot=[XYtot; XY{k}; NaN NaN];
    end
    difer(length(XYtot)-sum(lens)-length(lens),[],[],NaN)
	
    % Remove extra nans
    [xdata,ydata]=removeExtraNanSeparators(XYtot(:,1),XYtot(:,2));
    
    XYtot=[xdata ydata];
    
    % Make sure the coordinates make sense (0 to 360)
    XYtot(:,1)=XYtot(:,1)-360*[XYtot(:,1)>360];
    XYtot(:,1)=XYtot(:,1)+360*[XYtot(:,1)<0];
    % Run penlift
    [xdata,ydata]=penlift(XYtot(:,1),XYtot(:,2),200);
    %[xdata,ydata]=closePolygonParts(xdata,ydata,'degrees');
    XY=[xdata ydata];
    
    % Save the combined coordinates
    fnpl=sprintf('%s%s.mat',ddir2,'globe');
    save(fnpl,'XY')
  end
end

% Collect output
varns={XY};
% Provide output where requested
varargout=varns(1:nargout);

