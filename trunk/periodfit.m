function varargout=periodfit(yofx,perdx,x,comp)
% [vard,resd,pred,contr,perdb]=PERIODFIT(yofx,perdx,x,comp)
% 
% Performs a direct inversion for the periodic (sin/cos) components of a
% given signal - including an offset. But see POLYFIT for scaling on x.
%
% INPUT:
%
% yofx    The data, i.e. the dependent variable
% perdx   The target period, in the units of the independent variable
% x       The independent variable (if not present, equally spaced)
% comp    'series' searches for all periods individually, keeps best and
%                  returns the variance reduction for all tried
%         'parallel' searches for all periods at the same time
%
% OUTPUT:
%
% vard    Variance reduction of all periods tried
% resd    Residual of the best-fitting period
% pred    Prediction based on the best-fitting period
% contr   The estimated expansion coefficients (sin/cos)
% perdb   The best-fitting period from the set
% perdx   The target period if you didn't remember it
%
% Last modified by fjsimons-at-alum.mit.edu, 2/17/2012

defval('comp','series')

% Single-column data
yofx=yofx(:);

% Equally spaced data assumed by default
defval('x',[1:length(yofx)]')

% Faux output if nothing else
perdb=perdx;

if length(perdx)==1
  % The "design" matrix also includes an offset
  argm=x/perdx;
  F=[ones(size(yofx)) sin(2*pi*argm) cos(2*pi*argm)];
  % Do the inversion 
  [contr,pred,resd,vard]=doinversion(F,yofx);
else
  switch comp
   case 'series'
    % Recursive algorithm
    vard=nan(size(perdx));
    % Do them all
    for index=1:length(perdx)
      vard(index)=periodfit(yofx,perdx(index),x);
    end
    % Redo the best
    [~,vi]=max(vard);
    % Redefine the output to give us the best period
    perdb=perdx(vi);
    [~,resd,pred,contr]=periodfit(yofx,perdb,x);
   case 'parallel'
    % The "design" matrix also includes an offset
    argm=[1./perdx'*x']';
    F=[ones(size(yofx)) sin(2*pi*argm) cos(2*pi*argm)];
    
    % Do the inversion 
    [contr,pred,resd,vard]=doinversion(F,yofx);
   otherwise
    error('Pick a valid option, refer to the help')    
  end
end

% Optional output
varns={vard,resd,pred,contr,perdb,perdx};
varargout=varns(1:nargout);

% Do the actual inversion
function [contr,pred,resd,vard]=doinversion(A,y)

% Let's find the contributions from these periods
contr=inv(A'*A)*A'*y;

% Plot the prediction from this model
pred=A*contr;

% Let's look at the residuals
resd=pred-y;

% Compare the variance of the residual to the variance of the data: if
% the fit does a good job, the variance of the residuals should be a
% small fraction of the variance of the data
vard=100*[1-var(resd)/var(y)];
