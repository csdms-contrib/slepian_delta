function varargout=gracedeg1(Pcenter,Rlevel)
% [potcoffs,cal_errors,thedates]=GRACE2PLMT(Pcenter,Rlevel,units,forcenew)
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
%    potcoffs       potential coefficients [nmonths x addmup(Ldata) x 4]
%                    these could also be in surface mass density
%    thedates       time stamps in Matlab time
%
% NOTE:
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
% Last modified by charig-at-email.arizona.edu, 03/19/2021

% Determine parameters and set defaults
defval('Pcenter','CSR')
defval('Rlevel','RL06')



% Where the original data files are kept
defval('ddir1',fullfile(getenv('IFILES'),'GRACE','Degree1'));
 % Where you would like to save the new .mat file
defval('ddir2',fullfile(getenv('IFILES'),'GRACE','Degree1'));


% The name of the original data files
if strcmp(Rlevel,'RL06')
    fnpl1=sprintf('%s/TN-13_GEOC_%s_%s.txt',ddir1,Pcenter,Rlevel);
    fnpl2=sprintf('%s/TN-13_GEOC_%s_%s.mat',ddir2,Pcenter,Rlevel);
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
% automatically make and use a new .mat file if you update the source file.
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
     fclose(fid);
             
     % What we actually want is a 3 dimensional cell array with the date as
     % the first dimension and a lmcosi as dimension 2 and 3 (just like
     % grace2plmt)

     % Use the epoch start and end dates to get the midpoint date 
     thedates = (datenum(C{8},'yyyymmdd') + datenum(C{9},'yyyymmdd'))./2;
     [b,m] = unique(thedates);
     thedates = thedates(m);

     % Old formula
     % for i=1:n/2; temp = [deg1(2*i-1,2:7); deg1(2*i,2:7)]; mydeg1(i,:,:) = temp; end;
     
    % Open and scan the file (data from all three centers is 10 columns)
    %fid = fopen(fnpl1);
    %C = textscan(fid,'%s%s%s%s%s%s%s%s%s%s');


    %keyboard
    % Only grab the lines for GRCOF2
    Carray = cat(3,C{:});
    I = strmatch('GRCOF2',Carray(:,1,1),'exact');
    Carray = squeeze(Carray(I,1,:));
    
    % Only want columns 2-7, and as format double
    Carray = Carray(:,2:7);
    lmcosi_month=cellfun(@str2num,Carray);
    % This should be addmup(Ldata)
    [m,n] = size(lmcosi_month);
     
     
     
else
     load(fnpl)
     disp(sprintf('%s loaded by GRACEDEG1',fnpl))
end








% Collect output
varns={};
varargout=varns(1:nargout);




