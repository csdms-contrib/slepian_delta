function varargout=readgldas(version,resolution,timeres,forcing,varlist,forcenew)
% [datanames,varlist,thedates,thevars]=READGLDAS(version,resolution,timeres,forcing,varlist,forcenew)
%
% This program reads in the GLDAS products from NASA GSFC which should 
% be netcdf format, and saves them as matlab mat files.  It functions like
% a version of NET2MAT specific for GLDAS products.
%
%
% INPUT:
% 
% version       What version of GLDAS you want. [default: GLDAS1]
% resolution    The spatial resolution of the product you want in degrees. 
%                 Choose either:
%                 1 for one degree resolution [default]
%                 0.25 for 0.25 degree resolution
% timeres       The time resolution of the data you want. Choose either:
%                 'M' for monthly [default]
%                 '3H' for 3-hourly
% forcing       The model forcing used. [default: 'NOAH']
% varlist       The variables you want from the GLDAS products.  The
%                 default is all of them.
% forcenew      Wether or not you want to force new generation of a save file
%                (1) or just use the one we already have (0) [default].
%
% OUTPUT:
%
% datanames     The file names of the GLDAS files that were read.
% varlist       The list of variables back to you. Useful if you did not know
%               this beforehand and used the default.
% thedates      The dates corresponding to the files that you read.  If you
%               requested monthly data, then this date will be the midpoint
%               of each month since we know a priori the data averages over
%               an entire month.
% thevars       All of the variables you requested will be placed in a cell 
%               array in the this output. Each indicie in this cell array
%               is itself a cell array for each variable, for each file.  For 
%               example if you read two netcdf files, requesting two 
%               variables, then thevars will be 1x2, and thevars{1} will 
%               also be a 1x2 cell array.
%
% NOTES:
%   The GLDAS files can be found at http://ldas.gsfc.nasa.gov/gldas/index.php
%   One way to access them is through the Data and Information Services Center at
%   http://disc.sci.gsfc.nasa.gov/hydrology/data-holdings
%
%   Default behaviour is to save the variable mat files in the same
%   directory as the netcdf files themselves.
%
%
% EXAMPLE: 
%
% See also: NETVARREAD, NET2MAT, NCREAD
%
% Last modified by charig-at-princeton.edu, 03/16/2016

% Set defaults
defval('version','GLDAS1')
defval('forcenew',0)
defval('resolution',1)
defval('timeres','M')
defval('forcing','NOAH')
defval('timespan',1)
defval('data',1)
if resolution==0.25
    resolution = '025';
elseif resolution==1
    resolution = '10';
else
    error('There was a problem with the spatial reolution you requested.')
end

% Form a directory out of the model we want
if strcmp(timeres,'M')
    defval('ddir1',fullfile(getenv('IFILES'),'EARTHMODELS',version,[forcing resolution '_' timeres],filesep));
elseif strcmp(timeres,'3H')
    defval('ddir1',fullfile(getenv('IFILES'),'EARTHMODELS',version,[forcing resolution 'SUBP_' timeres],filesep));
else
    error('There was a problem with the time resolution you requested.');
end
    
% Check if this directory exists
if exist(ddir1,'dir')~=7
    error('The directory for the model you requested was not found.')
end

% Get the variable list if we do not have one
defval('varlist','netvarread(ls2cell(fullfile(ddir1,''GLDAS*001.nc''),1))');
if isstr(varlist)
    [varlist] = eval(varlist);
end

% Check if there are already saved .mat files there, corresponding to your
% varlist you gave
havefile=0;
for i=1:length(varlist)
    if exist([ddir1 varlist{i} '.mat'],'file'); % cannot use a cell array in the function exist
       havefile = havefile + 1;
    end
end
% Also need a dates file
if exist([ddir1 'thedates.mat'],'file'); havefile = havefile + 1; end

% Now read in the source files in this directory 
datanames=ls2cell(fullfile(ddir1,'GLDAS*001.nc'));

%%%
% The reading
%%%

% Do we have all the files we asked for? Load it.  Otherwise, or if we force it, make
% a new one (e.g. you added extra months to the database).
if havefile==length(varlist)+1 && forcenew==0
    for i=1:length(varlist)
       load([ddir1 varlist{i} '.mat'])
       disp(sprintf('%s loaded by READGLDAS',varlist{i}))
       
       % Now take all these separate variables and combine them into one
       % called thevars for output
       eval(['thevars{i} =' varlist{i} ';'])
    end
    load([ddir1 'thedates.mat'])
else
   % Make new files!
        
   % Initialize
   ndates = length(datanames);
   thedates = zeros(1,ndates);

   % Get "thedates" from the filenames
   if strcmp(timeres,'M')
       [thedates] = filename2date(datanames,{[17:22] 'yyyymm'});
   end

   % For GLDAS monthly data we know that they average an entire month.  So add
   % a half month to each date.
   if strcmp(timeres,'M')
       thedatesmid = thedates + eomday(indeks(datevec(thedates),':,1'),indeks(datevec(thedates),':,2'))'/2;
       thedates = thedatesmid;
   end

   % Read in all the files and same as mat files
   [datanames,thevars] = net2mat(ddir1,varlist,ddir1);
   
   % Create an extra save file which has thedates for the files we just read
   save([ddir1 'thedates.mat'],'thedates')

end

% Collect output
varns={datanames,varlist,thedates,thevars};
varargout=varns(1:nargout);









