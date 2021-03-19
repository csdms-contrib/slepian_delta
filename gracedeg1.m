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



d1 = dir(fnpl1);
d2 = dir(fnpl2);

% If this file already exists, and it is more new than the source file, 
% then load it.  Otherwise make a new one and return that one. This should
% automatically make and use a new .mat file if you update the source file.
if ~exist(fnpl2,'file')==2 || datenum(d1.date) > datenum(d2.date)
     % Here create a new save file to use.
     % This is reversed from our typical if, exist, load convention because
     % this way if the file does not exist it will skip the second if
     % condition (which would error on the datenum call) and go straight here
     
     
     
else
     load(fnpl)
     disp(sprintf('%s loaded by GRACEDEG1',fnpl))
end








% Collect output
varns={};
varargout=varns(1:nargout);




