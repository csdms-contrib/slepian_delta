function varargout=plm2avg(lmcosi,dom) 
% [Int,A,miniK,XY]=PLM2AVG(lmcosi,dom)
%
% Computes the integral and average value of a spherical harmonic 
% function (lmcosi) within a specific area (dom).  This is done by
% creating an integration vector and multiplying by the unwrapped
% coefficients of the field (from lmcosi). 
%
% INPUT:
%
% lmcosi    Standard-type real spherical harmonic expansion coefficients
% 
% dom       A string name with an approved region such as 'africa', OR
% XY        its coordinates (such as supplied from 'africa' etc)
%
% OUTPUT:
%
% Int       Integral of the function within the region
% A         Average value of the function within the region      
% miniK     The integration vector (also the first row/column from Klmlmp in KERNELC)
% XY        The coordinates of the region we integrated over
%
%
% EXAMPLE:
%
% plm2avg('demo1')  Integrate the output from GEOBOXCAP over a region and
%                   see that you get something close to the region's area
% plm2avg('demo2')  Compare this integration to summing the field which has
%                   been expanded onto an equal area Fibonacci grid
%
% SEE ALSO:
%   KERNELC, SPHAREA, FIBONACCI_GRID
%  
% Last modified by charig-at-princeton.edu, 06/27/2012

defval('dom','greenland')
defval('lmcosi',[0 0 0 0])

if ~isstr(lmcosi)

% Highest degree (bandwidth of the expansion)
Lmax=lmcosi(end,1);

% Get the coordinates
if isstr(dom)
  % Specify the region by name
  defval('pars',10);
  eval(sprintf('XY=%s(%i);',dom,pars));
else
  % Specify the region by the coordinates of the bounding curve
  XY=dom;
end
      
% Calculate northernmost and southernmost colatitude
thN=90-max(XY(:,2)); thN=thN*pi/180;
thS=90-min(XY(:,2)); thS=thS*pi/180;

% Introduce and dimensionalize variables and arrays
[dems,dels,mz,lmc,mzi,mzo,bigm,bigl,rinm,ronm]=addmon(Lmax);
dimK=(Lmax+1)^2; 
lenm=length(dems);

% Calculate some Gauss-Legendre points
intv=cos([thS thN]);
% Degree of Gauss-Legendre integration
ngl=200; 
nGL=min(ngl,size(XY,1)/2);
% These are going to be the required colatitudes - forget XY
[w,x,N]=gausslegendrecof(nGL,[],intv);
disp(sprintf('%i Gauss-Legendre points and weights calculated',N))
      
% Some arrays needed in a bit
% Note that dimK==sum(dubs)
dubs=repmat(2,lenm,1); dubs(mz)=1; comdi=[];

% Calculate all the Legendre functions themselves
Xlm=repmat(NaN,length(x),lenm);
ind=0;
for l=0:Lmax  
   Xlm(:,ind+1:ind+l+1)=[legendre(l,x(:)','sch')*sqrt(2*l+1)]';
   ind=ind+l+1;
end
%disp('Legendre functions calculated')

% No need for this to be longer, or to loop over anything
% Normally it would help take out the first dimension of redundancy between
% XlmXlmp and Klmlmp
comdi=dubs';
      
% In our ordering, the -1 precedes 1 and stands for the cosine term 
% This creates an indexing vector into Xlm that enlarges it from [0 01 012]
% to [0 0-11 0-11-22]
% Since we only have Xlm and not XlmXlm, we stop with coss.  In Kernelc,
% this was again expanded into bigo to account for this redundancy in two
% dimensions
comdex=[1:lenm];
coss=gamini(comdex,comdi);
bigo=coss;
      
% Get the longitudinal integration info for the domain
defval('Nk',10);
% Now we may have multiple pairs
phint=dphregion(acos(x)*180/pi,Nk,dom);
phint=phint*pi/180;

% No need to initialize miniK since we can make it all at once
ondex=0;
% Calculate the longitudinal integrals for each l,m combination
for lm1dex=1:dimK
   m1=bigm(lm1dex);
   ondex=ondex+1;
   if m1<=0
      I(:,ondex)=coscos(acos(x),m1,0,phint);
   elseif m1>0
      I(:,ondex)=sincos(acos(x),m1,0,phint);
   end
end
% coscos was reused here even though one of the m values is 0 because it is
% useful for doing odd geographics in matrix form (i.e. if the continent
% looks like a circle with hole in it)

% Make miniK, really the first row/column of Klmlmp from KERNELC
miniK = w(:)'*(Xlm(:,bigo).*I);
	  
% To make this exactly equivalent to Tony's \ylm, i.e. undo what we
% did above here, taking the output of YLM and multiplying
miniK=miniK/4/pi;

% Compare with KERNELC
%if xver==1
  %KK=rindeks(kernelc(Lmax,dom),1);
%end
     
% This then makes miniK(1) the fractional area on the sphere
% Check by comparing to output from spharea, if you want
% Note: SPHAREA defaults to 17 abcissas and weights, while this code uses 101.
% So differences to be expected when continents are squiggly.
% Check the first term which should equal the area on the unit sphere
%A1=spharea(XY); A2=areaint(XY(:,2),XY(:,1));
%disp(sprintf('Area check...  PLM2AVG A: %6.7f ; SPHAREA A: %6.7f ; AREAINT A: %6.7f',...
%	miniK(1),A1,A2))

% Now take the vector miniK and multiply with the vector for the function
% you want (lmcosi) and you will get the integral of that field within the
% region of interest
thecofs=lmcosi(:,3:4);
theINT=miniK*thecofs(mzo);

% To get the average value of the function in the region, divide by the
% area, which is likely most accurate in miniK (i.e. more accurate than SPHAREA)
A=theINT/miniK(1);

% Provide output where requested
varns={theINT,A,miniK,XY};
varargout=varns(1:nargout);

elseif strcmp(lmcosi,'demo1')
  % Integrate the output from GEOBOXCAP which puts 1 in the region and 0
  % elsewhere.  This should integrate to the area of the region.  In
  % practice, since the SH expansion of the GEOBOXCAP mask has rings, this
  % comparison has a good amount of error.
  % This is FRACTIONAL area on the unit sphere.
  dom='australia';
  L=40;
  degres=[];
  [Bl,dels,r,lon,lat,lmcosi]=geoboxcap(L,dom,[],degres,1);
  [Int,A,miniK,XY]=plm2avg(lmcosi,dom);
  maps=plotplm(lmcosi,[],[],4,degres); colorbar
  % Really dumb thing 
  %sum(maps(:))*2*pi/size(maps,2)*pi/size(maps,1)/4/pi 
  disp(' '); A1=spharea(XY);
  disp('Integration check... This should equal the area of the region with error mostly from GEOBOXCAP \n');
  disp(sprintf('PLM2AVG Int: %6.7f ; SPHAREA A: %6.7f ; ERROR: %6.2f%%',...
      Int,spharea(XY),(spharea(XY)-Int)/spharea(XY)*100))

  % Provide output where requested
  varns={Int,A,miniK,XY};
  varargout=varns(1:nargout);
elseif strcmp(lmcosi,'demo2')
  % Test the integration against expanding on a Fibonacci grid and summing
  dom='australia';
  Lmax=20;
  pars=10;
  disp(['Testing integration of EGM2008 (Lmax=' num2str(Lmax) ') over ' dom]);
  % Get coordinates
  eval(sprintf('XY=%s(%i);',dom,pars));
  % Make sure the coordinates make sense
  XY(:,1)=XY(:,1)-360*[XY(:,1)>360];
  % Use the first Lmax data from EGM2008 as a field of interest
  v=fralmanac('EGM2008_ZeroTide','SHM');
  % Note that gravity does not start at zero
  % Geoid = 3 % Free-air gravity anomaly = 2
  v=plm2pot(v(1:addmup(Lmax)-addmup(v(1)-1),:),[],[],[],3);
  % Create a Fibonacci grid
  [lonF,latF] = Fibonacci_grid(30000);
  % Expand EGM2008 on the Fib Grid
  [rF,lon,lat,Plm] = plm2xyz(v(:,1:4),latF,lonF);
  % Now decide if we're inside or outside of the region, and set outside
  % points to zero 
  rF(~inpolygon(lonF,latF,XY(:,1),XY(:,2)))=0;
  % Now compute the integration which can be 
  % discretized due to the equal area Fibonacci grid
  IntF = sum(rF)/length(rF);
  % Do the same integration with plm2avg
  [Int,A]=plm2avg(v(:,1:4),dom);
  % Now check the averages.
  indx = find(rF);
  AF=sum(rF(indx))/length(indx);
  
  disp(' ');
  disp('Integration check... PLM2AVG should equal the Fib Grid for good resolution');
  disp(sprintf('PLM2AVG Int: %6.7f ; FIB GRID Int: %6.7f ; ERROR: %6.2f%%',...
      Int,IntF,(Int-IntF)/Int*100))
  disp(' ');
  disp('Avg value check... PLM2AVG should equal the Fib Grid for good resolution');
  disp(sprintf('PLM2AVG Avg: %6.7f ; FIB GRID Avg: %6.7f ; ERROR: %6.2f%%',...
      A,AF,(A-AF)/A*100))
end

