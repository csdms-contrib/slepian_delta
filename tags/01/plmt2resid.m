function varargout=plmt2resid(plmt,thedates,fitwhat,givenerrors)
% [ESTresid,thedates,ESTsignal,rmset,rmsst,varet,varst]
%             =PLMT2RESID(plmt,thedates,fitwhat,givenerrors)
%
% Takes a time series of spherical harmonic coefficients and fits a desired
% combination of functions (e.g. secular, annual, semiannual, etc.)
%
% INPUT:
%
% plmt        The time series of spherical harmonic coefficients.  This
%               should be a three dimensional matrix (not a cell array), 
%               where the first dimension is time, and dimensions 2 and 3 
%               are as in the traditional lmcosi matrix.
% thedates    An array of dates corresponding to the plm thedateseries.  These
%               should be in Matlab's date format. (see DATENUM)
% fitwhat     The functions that you would like to fit to the time series
%               data.  The format of this array is as follows:
%               [mean secular periodic1 periodic2 periodic3 etc...] where
%               - mean is either 0/1 to turn off/on this function
%               - secular is either 0/1 to turn off/on this function
%               - periodic1 is the period in days of a function (i.e. 365.0)
%              Any # of desired periodic functions can be included. 
% givenerrors  These are given errors, if you have them.  In this case a
%                weighted inversion is performed.  givenerrors should be the
%                same dimensions of plmt.
%            
% OUTPUT:
%
% ESTresid    Residual time series for each pair of cos/sin coefficients
%             [nmonths x Ldata x 4] with a hyper-lmcosi format
%             ESTresid(:,:,1) the degree
%             ESTresid(:,:,2) the order
%             ESTresid(:,:,3) the cosine coefficient
%             ESTresid(:,:,4) the sine coefficient
% thedates    The dates that belong to these time series
% ESTsignal   The least-squares fitted function for each cos/sin
%             coefficient evaluated at those months in the same format 
% rmset       The root mean square of the estimated residual coefficients
% rmsst       The root mean square of the estimated signal coefficients
% varet       The variance of the estimated residual coefficients
% varst       The variance of the estimated signal coefficients
%
% SEE ALSO:
%
% Last modified by charig-at-princeton.edu 5/24/2011
% Last modified by fjsimons-at-alum.mit.edu 5/26/2011

defval('xver',0)  
defval('plmt',fullfile(getenv('IFILES'),'GRACE','CSR_alldata.mat'))
defval('plmt',fullfile(getenv('IFILES'),'GRACE','GFZ_alldata.mat'))
warning('off','MATLAB:divideByZero');

if isstr(plmt)
  % Load the file specified
  a=load(plmt);
  % These are the potential coefficients
  plmt=a.potcoffs(:,:,1:4);
  % These would be with the "formal" errors
  % givenerrros=a.potcoffs(:,:,[1 2 5:6]);
  % Much better to use the "calibrated" errors
  givenerrors=a.cal_errors(:,:,1:4);
  % Now get the thedates for each of those months
  thedates=a.thedates;
  % Get rid of the rest
  clear a
end

% Initialize/Preallocate
defval('givenerrors',ones(size(plmt)));
defval('fitwhat',[1 1 181.0 365.0])

% Round to 1 January of the first available and next last available year 
bounds=[datevec(thedates(1)); datevec(thedates(end))];
x0 = datenum(bounds(1)  ,1,1);
x1 = datenum(bounds(2)+1,1,1);
% The frequencies being fitted in [1/days]
omega = 1./[fitwhat(3:end)];
[i,j,k] = size(plmt);
% Initialize the residuals
ESTresid =zeros(size(plmt));
% Initialize the evaluated fitted function set
ESTsignal=zeros(size(plmt));
% Put in the degrees and the orders
ESTresid(:,:,1:2) =plmt(:,:,1:2);
ESTsignal(:,:,1:2)=plmt(:,:,1:2);
% Figure out the orders and degrees of this setup
Ldata=addmup(j,'r');
[dems,dels]=addmon(Ldata);
difer(squeeze(plmt(1,:,2))-dems(:)',[],[],NaN)
difer(squeeze(plmt(1,:,1))-dels(:)',[],[],NaN)
% How many data?
nmonths=length(thedates);
% How many periodic components?
lomega=length(omega);

% Assemble a G matrix for fitting
G = [];
% Mean term
if fitwhat(1) ~= 0
  G = [G ones(size(thedates'))];
end
% Secular term
if fitwhat(2) ~= 0
  G = [G (thedates-x0)'];
end
% Periodic terms
if ~isempty(omega)
  th= repmat(omega,nmonths,1)*2*pi.*repmat((thedates-x0)',1,lomega);
  G = [G cos(th) sin(th)];
end

% Initilization Complete

% Loop over coefficients starting from the (2,0) terms and do fitting
for index=4:j
    % If we have a priori error information, create a weighting matrix, and
    % change the G and d matrices to reflect this.  Since each coefficient
    % has its own weighting, we have to invert them separately.
    W3 = diag([1./squeeze(givenerrors(:,index,3))]);
    % Isolate the data over time
    d = squeeze(plmt(:,index,3:4));
    G3 = W3*G;
    if dems(index)==0
      % The sine coefficient and errors are 0. So prevent an error.
      G4 =  G;
      dw = [W3*d(:,1) d(:,2)];
    else
      W4 = diag([1./squeeze(givenerrors(:,index,4))]);
      G4 =  W4*G;
      dw = [W3*d(:,1) W4*d(:,2)];
    end
    
    % Do the fitting by minimizing least squares which is exactly right
    mL2 = [(G3'*G3)\(G3'*dw(:,1)) (G4'*G4)\(G4'*dw(:,2))];
    % mL2 is m by 2, where col 1 is for the C terms, col 2 for S terms
    
    % Use the model parameters to make the periodic amplitude in time
    startP = length(mL2) - 2*lomega + 1;
    Camp = [mL2(startP:(startP+lomega-1),1) mL2((startP+lomega):end,1)];
    Camp = sqrt(Camp(:,1).^2 + Camp(:,2).^2);
    
    Samp = [mL2(startP:(startP+lomega-1),2) mL2((startP+lomega):end,2)];
    Samp = sqrt(Samp(:,1).^2 + Samp(:,2).^2);
    
    % Use the model parameters to make the periodic phase in time
    Cphase = [mL2(startP:(startP+lomega-1),1) mL2((startP+lomega):end,1)];
    Cphase = atan2(Cphase(:,1),Cphase(:,2));

    Sphase = [mL2(startP:(startP+lomega-1),2) mL2((startP+lomega):end,2)];
    Sphase = atan2(Sphase(:,1),Sphase(:,2));
    
    % Assemble the estimated signal function, evaluated at 'thedates'
    fitfnC = 0;
    fitfnS = 0;
    
    % Start adding things
    if fitwhat(1) ~= 0 % We have a mean term
       fitfnC = fitfnC + mL2(1,1);
       fitfnS = fitfnS + mL2(1,2);
       if fitwhat(2) ~= 0 % And we have a secular term
           fitfnC = fitfnC + mL2(2,1)*(thedates-x0);
           fitfnS = fitfnS + mL2(2,2)*(thedates-x0);
       end
    else % We have no mean term, just perhaps a secular term
       if fitwhat(2) ~= 0
           fitfnC = fitfnC + mL2(1,1)*(thedates-x0);
           fitfnS = fitfnS + mL2(1,2)*(thedates-x0);
       end
    end
    

    % Add the sum over all the periodic components periodics
    fitfnC = fitfnC + ...
             sum(repmat(Camp,1,nmonths).*sin(th'+repmat(Cphase,1,nmonths)),1);
    fitfnS = fitfnS + ...
             sum(repmat(Samp,1,nmonths).*sin(th'+repmat(Sphase,1,nmonths)),1);
    
    % This is completely equivalent with the above for the purpose of residual
    % determination
    if xver==1
      % sum(repmat(mL2(startP:(startP+lomega-1),1),1,nmonths).*cos(th'),1)+...
      %       sum(repmat(mL2((startP+lomega):end,1),1,nmonths).*sin(th'),1)
      % sum(repmat(mL2(startP:(startP+lomega-1),2),1,nmonths).*cos(th'),1)+...
      %       sum(repmat(mL2((startP+lomega):end,2),1,nmonths).*sin(th'),1)
      clf
      subplot(211)
      plot(thedates,d(:,1),'+')
      hold on
      plot(thedates,fitfnC,'r');
      axis tight
      datetick('x',28)
      title(sprintf('cosine component at l = %i and m = %i',...
                    dels(index),dems(index)))
      
      subplot(212)
      plot(thedates,d(:,2),'+')
      hold
      plot(thedates,fitfnS,'r');
      axis tight
      datetick('x',28)
      title(sprintf('sine component at l = %i and m = %i',...
                    dels(index),dems(index)))
      pause(2)
    end
    
    % Compute the residual time series for this pair of cosine and sine
    % coefficients 
    resid = d - [fitfnC' fitfnS'];
    
    % Here's the definition of the residual at lm vs time
    ESTresid(:,index,3:4) = resid;
    % Here's the definition of the signal at lm vs time
    ESTsignal(:,index,3:4) = [fitfnC' fitfnS'];
end

% Compute the RMS of the estimated signal
rmsst = squeeze(ESTsignal(1,:,1:2));
rmsst(:,3:4) = squeeze( sqrt( mean( ESTsignal(:,:,3:4).^2,1 ) ) );

% Compute the RMS of the estimated residual
rmset = squeeze(ESTresid(1,:,1:2));
rmset(:,3:4) = squeeze( sqrt( mean( ESTresid(:,:,3:4).^2,1 ) ) );

% Compute the variance of the estimated signal
varst = squeeze(ESTsignal(1,:,1:2));
varst(:,3:4) = squeeze( var( ESTsignal(:,:,3:4),1 ) );

% Compute the variance of the estimated signal
varet = squeeze(ESTresid(1,:,1:2));
varet(:,3:4) = squeeze( var( ESTresid(:,:,3:4),1 ) );

% Collect output
varns={ESTresid,thedates,ESTsignal,rmset,rmsst,varet,varst};
varargout=varns(1:nargout);
