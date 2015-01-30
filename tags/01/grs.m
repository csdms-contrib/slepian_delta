function varargout=grs(GM,rf,a,omega,nzo) 
%========================= SYNTAX =================================
%
% [barC2n,geqt,gpol,U0,m,ecc,eccp,b,E,c]=
%                         grs(GM,rf,a,omega,nzo) 
%
%========================= PURPOSE ================================
% 
% Calculate geopotential constants in a reference earth model.
%
%========================= INPUTS =================================
%
% GM         Gravitational constant times mass reference   [scalar]
% rf         Inverse flattening (1/f)                      [scalar]
% a          Semi-major axis                               [scalar]
% omega      Angular velocity                              [scalar]
% nzo        Number of zonal harmonics (2,4,... 2*nzo)     [scalar]
%
%========================= OUTPUTS ================================
%
% barC2n     Normalized even zonal harmonics of            [matrix]
%            the corresponding Somigliana-Pizzetti normal field.          
%            barC2n(:,1): normalized zonal harmonics
%            barC2n(:,2): degree of the zonal harmonic [2 4 ... 2*nzo]
% geqt       Normal gravity at the equator                 [scalar]
% gpol       Normal gravity at the pole                    [scalar]
% U0         Normal potential at the ellipsoid             [scalar]
% m          omega^2*a^2*b/(GM)                            [scalar]
% ecc        First eccentricity                            [scalar]
% eccp       Second eccentricity                           [scalar]
% b          Semi-minor axis                               [scalar]  
% E          Linear eccentricity                           [scalar]
% c          Polar radius of curvature                     [scalar]
%
%========================= EXAMPLES ===============================
%
% [1] For the WGS84:
%       GM=0.3986004418e15;
%       a=6378137;
%       omega=7292115e-11;
%       rf=298.257223563;
%       [barC2n,geqt,gpol,U0,m,ecc,eccp,b,E,c]=grs(GM,rf,a,omega) 
% Compare the output to page in ref. [2].
%       
%========================= NOTES ==================================
%
% [1] All the following page numbers and equation numbers refer to the
%       book 'Physical Geodesy' by Hofmann-wellenhof and Moritz + 2006
%
%========================= REFERENCES =============================
%
% [1] Moritz (1984)
% [2] Hofmann-Wellenhof and Moritz (2006)
% [3] Heiskanen and Moritz (1964)
% [4] EGM2008: hsynth_WGS84.f
% 
%====== Last modified by dongwang-at-princeton.edu, 09/02/2009 ====
% Last modified by fjsimons-at-alum.mit.edu, 02/22/2012

% INPUT and OUTPUT error check.
error(nargchk(0,5,nargin,'struct'));
error(nargoutchk(0,10,nargout,'struct'));

% Default values for the WGS-84 ellipsoid
defval('GM',fralmanac('GM_wgs84'))
defval('rf',fralmanac('rf_wgs84'))
defval('a',fralmanac('a_wgs84'))
defval('omega',fralmanac('omega_wgs84'))
defval('nzo',10)

% Flattening
f=1/rf;

% First eccentricity;
try
  ecc=flat2ecc(1/rf);
catch
  ecc=sqrt(2/rf-rf.^-2);
end

% First eccentricity squared
ecc2=ecc^2;

% Second eccentricity
% p. 71, Eqn.(2-138);
eccp=ecc/(1-1/rf);
% Second eccentricity squared
eccp2=eccp^2;

% Semi-minor axis
b=a*(1-1/rf);

% m
% p. 70, Eqn.(2-137);
m=omega^2*a^2*b/GM;

% q_0 and q_0p
% p. 67, Eqn.(2-113)
% In the book, q_0 is q
q_0=1/2*((1+3/eccp2)*atan(eccp)-3/eccp);
% and q_0p is q_0
q_0p=3*(1+1/eccp2)*(1-1/eccp*atan(eccp))-1;

% J_2
% p. 75-76, Eqn.(2-165), Eqn.(2-166) and Eqn.(2-172)
j_2=ecc2/3*(1-2/15*m*eccp/q_0);

% j_2n
% p. 76, Eqn.(2-170) and Eqn.(2-172)
j_2n=J2N(nzo,ecc2,j_2);

% Normalized C2n0 terms.
% p. 60, Eqn.(2-80)
barC2n=zeros(nzo,2);
barC2n(:,2)=2*(1:nzo)';
barC2n(:,1)=-j_2n(:)./sqrt(barC2n(:,2)*2+1);

% Normal gravity at the equator.
% p. 71, Eqn.(2-141);
geqt=GM/a/b*(1-m-m/6*eccp*q_0p/q_0);

% Normal gravity at the pole.
% p. 71, Eqn.(2-142)
gpol=GM/a^2*(1+m/3*eccp*q_0p/q_0);

% Linear eccentricity
% p. 66, Eqn.(2-101)
E=sqrt(a^2-b^2);

% Normal potential at the ellipsoid
% p. 68, Eqn.(2-123)
U0=GM/E*atan(eccp)+1/3*omega^2*a^2;

% Polar radius of curvature 
% p. 73, eqn.(2-150)
c=a^2/b;

% Return things
vars={barC2n,geqt,gpol,U0,m,ecc,eccp,b,E,c};
varargout=vars(1:nargout);

function j_2n=J2N(nzo,E2,J2)
%========================= SYNTAX =================================
%
% j_2n=J2N(nzo,E2,J2)
%
%========================= PURPOSE ================================
% 
% Calculate the J_2N term (NOT normalized)
%
%========================= INPUTS =================================
%
% nzo         Index. 2*nzo is the even zonal harmonic.       [scalar]
% E2         First eccentricity squared                    [scalar]
% J2         J2 term.
%
%========================= OUTPUTS ================================
%
% J2         J_2 term.                                     [scalar]
%
%========================= NOTES ==================================
%
% [1] Here, J2 and j_2n are NOT normalized terms. The relationship
%         between normalized j_2n and un-normalized j_2n is given
%         by :
%            \bar{j_2n}=1/sqrt(2*2n+1)*j_2n
%
%========================= REFERENCES =============================
%
% [1] Heiskanen and Moritz (1964)
% [2] Hofmann-Wellenhof and Moritz (2006)
%
%====== Last modified by dongwang-at-princeton.edu, 09/02/2009 ====

% Make the index vector;
N=1:nzo;

% p. 76, Eqn.(2-170) and Eqn.(2-172)
j_2n=(-1).^(N+1).*3.*E2.^N./(2.*N+1)./(2.*N+3).*(1-N+5.*N.*J2./E2);

% Change to column vector
j_2n=j_2n(:);
