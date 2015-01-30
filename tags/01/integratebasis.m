function varargout=integratebasis(CC,TH,J)
% [eigfunINT]=INTEGRATEBASIS(CC,TH)
%
% Accepts a Slepian basis and integrates the functions within a region.
% The integrals reported are in fractional sphere area.  You may want to
% convert this to real sphere area via x(4*pi*radius^2).
%
% INPUT:
%
% CC      The eigenfunctions.  This can be either the G you get from 
%          GLMALPHA or a cell array of the functions as in LOCALIZATION. 
% TH      The region.  Can be formatted multiple ways as in GLMALPHA
% J       How many functions you want to do this for [DEFAULT: all]
%
% OUTPUT:
%
% eigfunINT    The integrals of the eigenfunctions over the region.  These
%               are in real sphere area, depending on your radius.
%
% SEE ALSO: PLM2AVG, PLM2AVGP
%
% Last modified by charig-at-princeton.edu, 06/27/2012

defval('TH','africa')  
defval('CC','[~,CC]=localization(15,TH);');

% Sort out what CC is
if isstr(CC)
  % Evaluate the specified expression
  eval(CC);
end
if isfloat(CC)
   % We must have a G matrix from glmalpha (already sorted)
   % Reorder them into a cell array
   [n,m] = size(CC);
   defval('J',m);
   L = sqrt(n) - 1;
   % This should be an integer
   if (floor(L)~=L);
      error('Something fishy about your L');
   end
   [~,~,~,lmcosi,~,~,~,~,~,ronm]=addmon(L);
   % Collect the eigenvector output into a format that PLM2XYZ knows how to interpret
   for j=1:J
      % Create the blanks
      cosi=lmcosi(:,3:4);
      % Stick in the coefficients of the 1st eigentaper
      cosi(ronm)=CC(:,j);
      % Construct the full matrix
      CC2{j} = [lmcosi(:,1:2) cosi]; 
   end
   CC = CC2;
elseif iscell(CC)
    % We are all good, just get the size.
    [n,m] = size(CC);
    defval('J',m);
    defval('L',addmup(size(CC{1},1),'r'));
else
    error('What format is CC?');
end

% Now sort out what TH is
if isstr(TH) % Geographic
    % No changes
    XY=TH;
elseif iscell(TH) % Geographic + buffer
    defval('buf',0);
    dom=TH{1}; buf=TH{2};
	defval('pars',10);
	% Run the named function to return the coordinates	
	XY=eval(sprintf('%s(%i,%f)',dom,pars,buf));
elseif isfloat(TH) && length(TH)==1 % Polar caps
    XY = [ [1:360]' repmat(90-TH,360,1) ]; 
end

% Initialization complete

%%%
% Do it
%%%

% See if we can run this calculation in parallel
s = matlabpool('size');
if s
  disp('Running INTEGRATEBASIS in parallel.')
  % plm2avgp is a parallel version, so no double parfor here
  for h=1:J
    [Int,A]=plm2avgp(CC{h},XY);
    eigfunINT(h) = Int;
  end
else
disp('Running INTEGRATEBASIS in series.  If unintended, open a matlabpool')

  % plm2avg is the regular version
  for h=1:J
    [Int,A]=plm2avg(CC{h},XY);
    eigfunINT(h) = Int;
  end
end

% Collect output
varns={eigfunINT};
varargout=varns(1:nargout); 



