function varargout=grace2plmt(Pcenter,Rlevel,units,forcenew)
% [potcoffs,cal_errors,thedates]=GRACE2PLMT(Pcenter,Rlevel,units,forcenew)
%
% This program reads in the Level-2 GRACE geoid products from either the CSR or
% GFZ data centers, does some processing, and saves them as a plmt matrix
% in a .mat file.  In particular, the coefficients are reordered to our
% prefered lmcosi format, they are referenced to the WGS84 ellipsoid, 
% the C2,0 coefficients are replaced with more accurate measurements from
% satellite laser ranging from Loomis et al, (2019), and the degree one coefficients are 
% substituted with those from Sun et al., (2016).  You have the option 
% of leaving them as geopotential 
% or converting them to surface mass density using the method of 
% Wahr et al. 1998, based on Love numbers (see PLM2POT).
%
% INPUT:
% 
% Pcenter     'CSR' data center at the Center for Space Research
%             'GFZ' data center at the GeoForschungsZentrum Potsdam
% Rlevel      The release level of the solution you want.  
%              Either 'RL04','RL05', or 'RL06'
% units       'POT' or 'SD' for whether you want geopotential or surface
%               mass density
% forcenew    Whether or not you want to force new generation of a save file
%              (1) or just use the one we already have (0) [default].
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
%   2/19/2021 Formal or calibrated uncertainties have not been reported 
%    since RL04 so we will discontinue the output of these. The output
%    variables will be altered so existing scripts will need to be adjusted
%    to take this into account. A corresponding change has been made in
%    GRACE2SLEPT
%
%	SLR data available from the GRACE Tellus website for RL06:
%	ftp://podaac-ftp.jpl.nasa.gov/allData/grace/L2/CSR/RL06
%	 
%
%   SLR data available from the GRACE Tellus website:
%   http://grace.jpl.nasa.gov/data/J2/ notably 
%   ftp://ftp.csr.utexas.edu/pub/slr/degree_2/C20_RL04_2010_12.txt
%  The header was removed and the file renamed for easy use.
%  Updated files keep getting posted in the same location.  
%
%  SLR data available from the GRACE Tellus website:
%   http://grace.jpl.nasa.gov/data/degree1/ notably 
%   ftp://podaac.jpl.nasa.gov/allData/tellus/L2/degree_1/
%  The header was removed and the file renamed for easy use.
%  Updated files keep getting posted in the same location.  
%
% EXAMPLE: to make a new save file when you have added more months
% [potcoffs,thedates]=grace2plmt('CSR','RL05','SD',1);
%
% Last modified by charig-at-email.arizona.edu, 02/19/2021
% Last modified by lashokkumar-at-arizona.edu, 11/09/2020
% Last modified by mlubeck-at-email.arizona.edu, 03/18/2019
% Last modified by fjsimons-at-alum.mit.edu, 05/17/2011

% Determine parameters and set defaults
defval('Pcenter','CSR')
defval('Rlevel','RL06')
defval('units','SD')
defval('forcenew',1)