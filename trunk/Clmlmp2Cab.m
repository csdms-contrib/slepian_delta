function varargout=Clmlmp2Cab(Clmlmp,TH,eigs)
% [Cab,Cab_scaled,G,V]=CLMLMP2CAB(Clmlmp,TH,eigs)
%
% Takes a spherical harmonic covariance matrix and calculates the
% covariance of a desired Slepian basis. (Cab = G*Clmlmp*G')  The Slepian
% basis will be the same bandwidth as the spectral covariance data
%
% INPUT:
%
% Clmlmp      Spherical harmonic covariance matrix, assumed to be ordered as
%              in GLMALPHA (not KERNELC)
% TH          The localization domain (a region string like 'greenland' 
%              OR [lon lat] an ordered list defining a closed curve [degrees]
%              OR a cell containing a region and a buffer such 
%               as {'greenland' 0.5})
%              
% eigs        Number of eigenfunctions you want in the basis
%               0 - up to the Shannon number [default]
%
% OUTPUT:
%
% Cab         Slepian basis covariance matrix
% Cab_scaled  Scaled Slepian basis covariance matrix
% G           The unitary matrix of localization coefficients (sorted)
% V           The eigenvalues (sorted)
%
% SEE ALSO:
%
% Last modified by charig-at-princeton.edu 6/26/2012

% Determine parameters and set defaults
defval('Lwindow',20)
defval('TH','greenland')
defval('pars',10)

% Top level directory
% For Chris
IFILES=getenv('IFILES');
% For FJS, who has a different $IFILES
%IFILES='/u/charig/Data/';

% Determine the bandwidth of our data.  Assume Clmlmp starts at l=0.
[n,m]=size(Clmlmp);
Ldata=addmoff(n,'r');
% J is how many Slepian functions we want to see
defval('eigs',round((Ldata+1)^2*spharea(TH)))

% Get a Slepian basis
[G,V,EL,EM,N,GM2AL,MTAP,IMTAP]=glmalpha(TH,Ldata,[],0,[],[],eigs);

% Sort by decreasing eigenvalue
[V,vi]=sort(V,'descend');
G=G(:,vi); if ~isnan(MTAP); MTAP=MTAP(vi); end
% If you don't do this, the eigenfunctions are ordered in the way
% that they correspond to single-orders back when, unrotated, they
% belonged to a polar cap, and the eigenvalues are sorted within
% these blocks. This is useful for, e.g. SPIE2009_1 a la SDSNEEUW.

% Slepian domain covariance
Cab = G'*Clmlmp*G;
% Scaled version
Cab_scaled=Cab./[diag(sqrt(Cab))*diag(sqrt(Cab))'];

% Collect output
varns={Cab,Cab_scaled,G,V};
varargout=varns(1:nargout);
