function varargout=gldas2TWSt(version,gldasres,timeres,forcing,masked,forcenew)
% [thedates,TWSlmcosi,TWS,lonlon,latlat]=GLDAS2TWST(version,gldasres,timeres,forcing,forcenew)
%
% This program reads in the GLDAS products from NASA GSFC which we have saved.
% It then takes the necessary variables and constructs a field for
% Terrestrial Water Storage for each month.  Then we remove areas where glaciers exist
% since these products are known to be inaccurate in these areas.  Finally
% we use xyz2plm to make spherical harmonics.
%
%
% INPUT:
% 
% version       What version of GLDAS you want. [default: GLDAS1]
% gldasres      The spatial resolution of the product you want in degrees. 
%                Choose either 1 or Currently only [default: 1]
% timeres       The time resolution of the data you want. Choose either:
%                'M' for monthly [default]
%                '3H' for 3-hourly
% forcing       The model forcing used. [default: 'NOAH']
% varlist       The variables you want from the GLDAS products.  The
%                default is all of them.
% masked        Do you want glacier areas masked? [default: 1 (yes)]
% forcenew      Wether or not you want to force new generation of a save file
%               (1) or just use the one we already have (0) [default].
%
% OUTPUT:
%
% thedates    The dates of the GLDAS files that were used
% TWSlmcosi   The lmcosi matrices for Terrestrial Water Storage for each 
%             month as a cell array
% TWS         Terrestrial Water Storage for each month as a spatial matrix
% lonlon      Longitude corresponding to the matrices in TWS
% latlat      Latitude corresponding to the matrices in TWS
%
%
% EXAMPLE: 
%
% See also: READGLDAS, GETGLACIERS, MASKFROMGMT, COARSENMASK, NET2MAT
%
% Last modified by charig-at-princeton.edu, 03/16/2016

% Set some defaults
defval('version','GLDAS1');
defval('timeres','M');
defval('forcing','NOAH');
defval('gldasres',1);
defval('c11cmn',[0.5 89.5 359.5 -89.5]);
defval('gldasc11cmn',[0.5 89.5 359.5 -59.5]);
defval('forcenew',0);
defval('masked',1);
% Any variables used in a parfor loop must be explicitly assigned during a
% function call.  i.e. using the ASSIGNIN function will not work
TWSpadding = [30 360];
L = 60;

% Some specific handling for the spatial resolution
if gldasres==0.25
    sres = '025';
elseif gldasres==1
    sres = '10';
else
    error('There was a problem with the spatial reolution you requested.')
end

% Where you would like to save the new .mat file
defval('ddir1',fullfile(getenv('IFILES'),'EARTHMODELS','GLDAS1'));
% Get the file name to save
if masked
fnpl=sprintf('%s/TWS-%s-%s-%s-%s-%s.mat',ddir1,version,forcing,timeres,sres,'masked');
else
fnpl=sprintf('%s/TWS-%s-%s-%s-%s-%s.mat',ddir1,version,forcing,timeres,sres);    
end
% If this expansion already exists, load it.  Otherwise, or if we force 
% it, make a new one (e.g. if you added extra months to the database).
if exist(fnpl,'file')==2 && forcenew==0
     load(fnpl)
     disp(sprintf('%s loaded by GLDAS2TWST',fnpl))
else

    % Get the GLDAS products to make TWS
    getvars = {'g0_lon_1' 'g0_lat_0' 'SoilMoist_GDS0_DBLY' 'Canopint_GDS0_SFC' 'SWE_GDS0_SFC'};
    [datanames,varlist,thedates,thevars] = readgldas(version,gldasres,timeres,forcing,getvars,forcenew);

    % Make Terrestrial Water Storage
    % These variables all have units of kg/m^2, and we want their sum
    for i=1:length(thedates)
        %totalsoil = sum(thevars{3}{i},3);
        TWS{i} = sum(thevars{3}{i},3) + thevars{4}{i} + thevars{5}{i};
    end

    % NOTE: Even though we mask this for glacier areas (areas where GLDAS is 
    % thought to be inaccurate), there are still problematic areas in Tibet.
    % Make sure to check this later

    if masked
      % Get Glacier/Antarctica/Greenland coordinates
      XYg = getglaciers('globe');
      XYgreen = greenland(10,0.5);
      XYant = antarctica(10,0.5,1);
      % Need to make Ant work
      [latF,lonF] = flatearthpoly(XYant(:,2),XYant(:,1),180);
      XYant = [lonF latF];
      % Combine these
      XYtot = [XYg; NaN NaN; XYant; NaN NaN; XYgreen];

      % Make the mask for these places
      [mymask,reshapedmask] = maskfromgmt(XYtot,num2str(gldasres/10),...
          [num2str(gldasc11cmn(1)) '/' num2str(gldasc11cmn(3)) '/' num2str(gldasc11cmn(4)) '/' num2str(gldasc11cmn(2))]);
      [mymask,lonlon,latlat] = coarsenmask(gldasres,gldasc11cmn,mymask,1);

      % Switch the sense of the mask, and make the lat/lon global for later
      mymask = (mymask-1)*-1; % Now glaciers are 0, not glaciers is 1
      coarselon = [c11cmn(1):gldasres:c11cmn(3)];
      coarselat = [c11cmn(2):-gldasres:c11cmn(4)];
      [lonlon,latlat]=meshgrid(coarselon,coarselat);
    else
      coarselon = [gldasc11cmn(1):gldasres:gldasc11cmn(3)];
      coarselat = [gldasc11cmn(2):-gldasres:gldasc11cmn(4)];
      mymask = ones(length(coarselat),length(coarselon)); 
      [lonlon,latlat]=meshgrid([c11cmn(1):gldasres:c11cmn(3)],[c11cmn(2):-gldasres:c11cmn(4)]);
    end

    % Multiply the mask to Total Water Storage in order to remove areas 
    % with permanent land ice.
    parfor i=1:length(thedates)
        TWS{i} = mymask.*circshift(rot90(TWS{i}),[0 180]);
%        TWS{i} = circshift(rot90(TWS{i}),[0 180]);
        TWS{i}(isnan(TWS{i}))=0;
        TWS{i} = [TWS{i}; zeros(TWSpadding(1),TWSpadding(2))];
        TWS{i} = [TWS{i} TWS{i}(:,1)];
        TWSlmcosi{i} = xyz2plm(TWS{i},L,[],latlat(:,1)');
        TWSlmcosi{i} = plm2rot(TWSlmcosi{i},-0.5);
    end
    % Frederik, is there a way to do a regular grid lmcosi inversion when the
    % data do not exactly end on 0 and 360?  This grid is basically global.
    
    % SAVE
    save(fnpl,'thedates','TWSlmcosi','TWS','lonlon','latlat');
    
end


% Collect output
varns={thedates,TWSlmcosi,TWS,lonlon,latlat};
% Provide output where requested
varargout=varns(1:nargout);








