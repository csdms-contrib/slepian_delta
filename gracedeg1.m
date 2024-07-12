function varargout=gracedeg1(Pcenter,Rlevel)
% [thedates,deg1data]=GRACEDEG1(Pcenter,Rlevel)
%
% This function reads and formats the degree-1 spherical harmonic
% correction from GRACE/GRACE-FO Technical Note 13.
%
% INPUT:
%
% Pcenter     'CSR' data center at the Center for Space Research
%             'GFZ' data center at the GeoForschungsZentrum Potsdam
%             'JPL' data center at the Jet Propulsion Laboratory
% Rlevel      The release level of the solution you want.
%              Either 'RL04','RL05', or 'RL06'
%
% OUTPUT:
%
% Returns these variables and saves them in a .mat file:
%    thedates       time stamps in Matlab time. This date is the midpoint
%                      of the GRACE month data span.
%    deg1data       potential coefficients for just degree l = 1
%                      [nmonths x 2 x 4]
%
% NOTES:
%
% The GRACE degree-1 correction data are distributed in Technical Note 13,
%   TN-13. There is one for each datacenter because the GRACE coefficients
%   from degree 2-60 are used in the construction of degree one terms. You
%   can find TN-13 in the docs folder in the GRACE data directory at NASA
%   PODAAC. https://podaac.jpl.nasa.gov/
%
% We check the file timestamps of both the TN-13 file and the Matlab .mat
%   file that we save/load. If the .mat file is newer, then we load it. If
%   the TN-13 document is newer, then we use it to create a new .mat file
%   and then load the .mat file.
%
% Last modified by charig-at-arizona.edu, 11/18/2021

% Determine parameters and set defaults
defval('Pcenter','CSR')
defval('Rlevel','RL06')



% Where the original data files are kept
defval('ddir1',fullfile(getenv('IFILES'),'GRACE','Degree1'));
% Where you would like to save the new .mat file
defval('ddir2',fullfile(getenv('IFILES'),'GRACE','Degree1'));


% The name of the original data files
if strcmp(Rlevel,'RL06')
    fnpl1=sprintf('%s/TN-13_GEOC_%s_%s.txt',ddir1,Pcenter,'RL0602');
    fnpl2=sprintf('%s/TN-13_GEOC_%s_%s.mat',ddir2,Pcenter,'RL0602');
elseif strcmp(Rlevel,'RL05')
    fnpl1=sprintf('%s/deg1_RL05_NH.txt',ddir1);
    fnpl2=sprintf('%s/deg1_RL05_NH.mat',ddir2);
elseif strcmp(Rlevel,'RL04')
    fnpl1=sprintf('%s/deg1_RL04_NH.txt',ddir1);
    fnpl2=sprintf('%s/deg1_RL04_NH.mat',ddir2);
else
    error('GRACEDEG1: Wonky release level requested')
end


% Get the file date information
d1 = dir(fnpl1);
d2 = dir(fnpl2);

% If this file already exists, and it is more new than the source file,
% then load it.  Otherwise make a new one and return that one. This should
% automatically make and use a new .mat file if you update the source (TN-13) file.
if ~(exist(fnpl2,'file')==2) || datenum(d1.date) > datenum(d2.date)
    % Here create a new save file to use.
    % This is reversed from our typical if, exist, load convention because
    % this way if the file does not exist it will skip the second if
    % condition (which would error on the datenum call) and go straight here
        
    fid = fopen(fnpl1);
    tline = fgetl(fid);
    
    % Keep reading lines until we get to the 'end of header' line
    while isempty(tline) || ~strcmp(tline(1:5),'end o')
        tline = fgetl(fid);
    end
    
    % Now do a formatted text scan on the data
    C = textscan(fid,'%s%d%d%f%f%f%f%s%s%s');
    % Use string format for the dates so we can easily send it to datenum
    fclose(fid);
    
    % What we actually want is a 3 dimensional cell array with the date as
    % the first dimension and a lmcosi as dimension 2 and 3 (just like
    % grace2plmt)
    
    % Use the epoch start and end dates to get the midpoint date
    thedates = (datenum(C{8},'yyyymmdd') + datenum(C{9},'yyyymmdd'))./2;
    [b,m] = unique(thedates);
    thedates = thedates(m);
    
    % Re cast the ints as doubles to get it all in one matrix
    deg1data = [double(C{2}) double(C{3}) C{4} C{5}];
    
    % Now reorder into a lmcosi format for each month
    for i=1:(length(C{1}))/2
        temp = [deg1data(2*i-1,:);
            deg1data(2*i,:)];
        deg1lmcosi(i,:,:) = temp;
    end
    
    deg1data = deg1lmcosi;
    
    % Create a save file
    save(fnpl2,'thedates','deg1data')
    
else
    load(fnpl2)
    fprintf('%s loading %s\n', upper(mfilename), fnpl2)
end





% Collect output
varns={thedates,deg1data};
varargout=varns(1:nargout);




