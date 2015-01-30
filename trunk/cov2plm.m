function varargout=cov2plm(Clmlmp,howmany)
% [lmcosiS]=COV2PLM(cov,howmany)
%
% This function accepts a covariance matrix for spherical harmonics and
% determines a set of data on the sphere (in SH coefficients) that has the
% properties of the covariance matrix you input.
%
% INPUT:
%
% Clmlmp        The covariance matrix of a set of spherical harmonics.
% howmany       The number of noise realizations you want [default: 1]
%            
% OUTPUT:
%
% lmcosiS     This will be a three dimensional matrix.  The first and 
%              second dimensions are or familiar matrix of spherical 
%              harmonic coefficients that represent a noise realization.
%              The third dimension is the number of noise realizations you 
%              wanted.  If you only wanted one set of data, then this
%              matrix will only have two dimensions, however you code will
%              still work because if you call lmcosiS(5,:,1) it is the same
%              as lmcosiS(5,:) for 2D matrices.
%
% SEE ALSO:
%
% EXAMPLE:
% lmcosi=cov2plm(Clmlmp,2);  plotplm(lmcosi(:,:,1)); figure; plotplm(lmcosi(:,:,2));
% 
% Last modified by charig-at-princeton.edu 06/26/2012

defval('Clmlmp','grace2plmt(''CSR'',''SD'')');
defval('xver',0);
defval('howmany',1);

if isstr(Clmlmp)
  % Evaluate the specified expression, and calculate a new covariance matrix
  [plmt,givenerrors,thedates] = eval(Clmlmp);
  [ESTresid]=plmt2resid(plmt,thedates,[1 1 365.0],givenerrors);
  [Clmlmp]=plmresid2cov(ESTresid);
end

% Decompose the covariance
T = cholcov(Clmlmp);
[n,m] = size(T);

% Check if this is right
if xver
    % Generate a lot of data that averages to the correct covariance 
    % (aside from random variation).
    SYNClmlmp = cov(randn(10000,n)*T);
    % Display some elements, and see if they are close to eachother
    Clmlmp(1:10,1:10)
    SYNClmlmp(1:10,1:10)
    pause
end

% How much data do we have?
Ldata = addmoff(m,'r');
% Get empty matrix and info for the data bandlimit
[~,~,~,lmcosidata,~,~,~,~,~,ronmdata]=addmon(Ldata);

% Make a synthetic noise realization
for i=1:howmany
   syntheticnoise = randn(1,n)*T;
   % Reorder the noise into lmcosi
   temp1=lmcosidata(:,3:4);
   temp1(ronmdata) = syntheticnoise(:);
   lmcosiS(:,:,i) = [lmcosidata(:,1:2) temp1];
end

% Reduce dimension if you only wanted one
lmcosiS = squeeze(lmcosiS);

% Collect output
varns={lmcosiS};
varargout=varns(1:nargout);
