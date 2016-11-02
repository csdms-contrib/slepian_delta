function varargout=makenormalmatrix(xdata,fitwhat)
% [N,mu]=MAKENORMALMATRIX(thedates,fitwhat,givenerrors)
%
% Takes a x-vector and information on what you plan to fit, and creates a
% normal matrix used to perform an inversion, e.g. least squares fitting.
%
% You can choose to fit either a mean, linear, quadratic, or cubic fuction
% as your "base" function to each Slepian coefficient, by using "fitwhat".
% In these cases of higher functions, they are only used if they pass an
% F-test for variance reduction.  For example, a cubic term is only used if
% it is significantly better than a quadratic function.
%
% INPUT:
%
% xdata       An array of xdata.
% fitwhat     The functions that you would like to fit to the time series
%               data.  The format of this array is as follows:
%               [order periodic1 periodic2 periodic3 etc...] where
%               - order is either 0/1/2/3 to fit UP TO either a
%                  mean/linear/quadratic/cubic function (if increasing
%                  the power of the function reduces variance enough)
%               - periodic1 is the period in the same units as your data 
%                  (e.g. 365.0 days for data in days)
%               Any # of desired periodic functions [days] can be included. 
%            
% OUTPUT:
%
% N           The normal matrix to be used in the inversion 
% mu          The scaling parameters used to adjust the xdata to make the
%               inversion more numerically stable. First is xdata mean and
%               second is xdata standard deviation, as in polyfit.
%
% SEE ALSO: POLYFIT
%
% Last modified by charig-at-princeton.edu  11/01/2016


% Initialize/Preallocate
defval('fitwhat',[3 365.0])
defval('xdata',[1:10]);

% How many data?
thismany=length(xdata);
xdata = xdata(:);
% We will do a scaling to improve the solution
mu1 = mean(xdata); % mean
mu2 = std(xdata); % standard deviation
% Make a new x-vector with this information
xprime = (xdata - mu1)/mu2;

mu = [mu1 mu2];

% The frequencies being fitted in [1/days]
omega = 1./[fitwhat(2:end)];
% Rescale these to our new xprime
omega = omega*mu2;

% How many periodic components?
lomega=length(omega);

%%%
% N matrix assembly
%%%

% The first section is basically a Vandermonde matrix, so we can make one
% flip it, and continue onwards.
N(:,1) = ones(length(xdata),1);
for j = 2:1:fitwhat(1)+1
   N(:,j) = xdata.*N(:,j-1);
end

% Periodic terms
if ~isempty(omega)
  % Angular frequency in radians/(rescaled day) of the periodic terms
  th_o= repmat(omega,thismany,1)*2*pi.*repmat((xprime),1,lomega);
  N = [N cos(th_o) sin(th_o)];
end 

% Collect the expanded output
varns={N,mu};

varargout=varns(1:nargout);
