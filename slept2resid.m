function varargout=slept2resid(slept,thedates,fitwhat,givenerrors,specialterms,CC,TH)
% [ESTsignal,ESTresid,ftests,extravalues,total,alphavarall,totalparams,
%    totalparamerrors,totalfit,functionintegrals,alphavar]
%       =SLEPT2RESID(slept,thedates,fitwhat,givenerrors,specialterms,CC,TH)
%
% Takes a time series of Slepian coefficients and fits a desired
% combination of functions (e.g. secular, annual, semiannual, etc.) to 
% each coefficient in order to separate "signal" from residual.
%
% You can choose to fit either a mean, linear, quadratic, or cubic fuction
% as your "base" function to each Slepian coefficient, by using "fitwhat".
% In these cases of higher functions, they are only used if they pass an
% F-test for variance reduction.  For example, a cubic term is only used if
% it is significantly better than a quadratic function.
%
% If you also provide the Slepian functions (CC) and region (TH) for this
% localization, then we assume that you want the fitting to the integrated
% functions. i.e. If you have surface density, then using the integrated
% Slepian function would give you total mass change in a region.  In
% addition, there will be a fitting of the combination (total) of the Slepian
% functions up to the Shannon number for this localization.
%
% INPUT:
%
% slept       The time series of Slepian coefficients.  This
%               should be a two dimensional matrix (not a cell array), 
%               where the first dimension is time, and the second dimension
%               are Slepian coefficients sorted by global eigenvalue. 
% thedates    An array of dates corresponding to the slept timeseries.  These
%               should be in Matlab's date format. (see DATENUM in units
%               of decimal days)
%             ***It is assumed that 'thedates' matches 'slept'.  If 'thedates'
%              is longer than 'slept' then it is assumed that these are
%              extra dates that you would like ESTsignal to be evaluated at.
%              This is useful if you want to interpolate to one of the
%              GRACE missing months that is important, such as Jan 2011.
% fitwhat     The functions that you would like to fit to the time series
%               data.  The format of this array is as follows:
%               [order periodic1 periodic2 periodic3 etc...] where
%               - order is either 0/1/2/3 to fit UP TO either a
%                  mean/linear/quadratic/cubic function (if increasing
%                  the power of the function reduces variance enough)
%               - periodic1 is the period in days of a function (i.e. 365.0)
%              Any # of desired periodic functions [days] can be included. 
% givenerrors  These are given errors, if you have them.  In this case a
%                weighted inversion is performed.  givenerrors should be the
%                same dimensions of slept.
% specialterms  A cell array such as {2 'periodic' 1460}.  At the moment
%                this is pretty specific to our needs, but could be 
%                expanded later.
% CC           A cell array of the localization Slepian functions
% TH           The region (proper string or XY coordinates) that you did the
%               localization on (so we can integrate)
%            
% OUTPUT:
%
% ESTsignal   The least-squares fitted function for each Slepian
%             coefficient evaluated at those months in the same format 
% ESTresid    Residual time series for each Slepian coefficients, ordered 
%              as they were given, presumably by eigenvalue
%              [nmonths x (Lwindow+1)^2] 
% ftests       An matrix, such as [0 1 1] for each Slepian coefficient, 
%                 on whether the fits you
%                 requested passed an F-test for significance.
% extravalues  These are the values of ESTsignal evaluated at your extra
%                 dates which were tacked onto 'thedates'
% total        The time series of the combined mass change from N Slepian
%                functions (i.e. the combined data points)
% alphavarall  The time averaged variance on each data point (from error 
%                propogation from each individual function).  These values
%                are the same for every point.
%
% totalparams   The parameters of the fits to the total.  This is a 4-by-j
%                 matrix where j are the fits you wanted.  Zeros fill out 
%                 the unused parameters.  For example, if
%                 you want just a 2nd order polynomial, you will get 
%                 [intercept intercept; slope slope; 0 quadratic; 0 0]
% totalparamerrors      The 95% confidence value on this slope.  This is computed
%               using the critical value for a t-distribution
%               At the moment, totalparams and totalparamerrors and just 
%               the values for a linear fit.  Sometime later maybe change 
%               this to be potentially a quadratic fit as well.
%
% totalfit       Datapoints for the best fit line, so you can plot it, and
%                 datapoints for the confidence intervals of this linear fit,
%                 so you can plot them.  NOTE: To plot these you should use
%                 totalfit(:,2)+totalfit(:,3) and totalfit(:,2)-totalfit(:,3)
%
% functionintegrals   This is a vector of the integrals of the Slepian
%                      functions, up to the Shannon number.  With this we
%                      can multiply by the data or ESTSignal to get the
%                      changes in total mass over time for each function.
%                      This also has the conversion from kg to Gt already
%                      in it, so don't do that again.
%
% alphavar    The variances on the individual alpha function time series (INTEGRALS).
%              This is calculated by making a covariance matrix from
%              ESTresid and then doing the matrix multiplication with the
%              eigenfunction integrals.
%
% SEE ALSO:
%
% Last modified by charig-at-princeton.edu  6/26/2012

defval('xver',0)  
%defval('specialterms',{2 'periodic' 1728.1});
defval('specialterms',{NaN});
defval('slept','grace2slept(''CSR'',''greenland'',0.5,60,[],[],[],[],''SD'')')
defval('extravalues',[]);

if isstr(slept)
  % Evaluate the specified expression
  [slept,~,thedates,TH,G,CC,V] = eval(slept);
end

% Initialize/Preallocate
defval('givenerrors',ones(size(slept)));
defval('fitwhat',[3 365.0])

% Handle the dates
if length(thedates)==size(slept,1)
    % Do nothing
    extradates = []; moredates=0;
elseif length(thedates) > size(slept,1)
    % There are extra dates, pull them off
    extradates = thedates((size(slept,1)+1):end);
    thedates = thedates(1:size(slept,1));
    moredates=1;
else
   error('Is thedates shorter than slept?') 
end

% How many data?
nmonths=length(thedates);
% We will do a scaling to improve the solution
mu1 = mean(thedates); % mean
mu2 = std(thedates); % standard deviation
% Make a new x-vector with this information
xprime = (thedates - mu1)/mu2;
extradatesprime = (extradates - mu1)/mu2;

% The frequencies being fitted in [1/days]
omega = 1./[fitwhat(2:end)];
% Rescale these to our new xprime
omega = omega*mu2;
[i,j] = size(slept);
% Initialize the residuals
ESTresid =zeros(size(slept));
% Initialize the evaluated fitted function set
ESTsignal=zeros(size(slept));
% Figure out the orders and degrees of this setup
% BUT orders and degrees have lost their meaning since slept should
%  be ordered by eigenvalue
Ldata=addmoff(j,'r');
[dems,dels]=addmout(Ldata);
% How many periodic components?
lomega=length(omega);

%%%
% G matrix assembly
%%%

% We will have the same number of G matrices as order of polynomial fit.
% These matrices are smallish, so make all 3 regardless of whether you want
% them all.
G1 = []; % For line fits
G2 = []; % For quadratic fits
G3 = []; % For cubic fits
% Mean term
if fitwhat(1) >= 0
  G1 = [G1 ones(size(xprime'))];
  G2 = [G2 ones(size(xprime'))];
  G3 = [G3 ones(size(xprime'))];
end
% Secular term
if fitwhat(1) >= 1
  G1 = [G1 (xprime)'];
  G2 = [G2 (xprime)'];
  G3 = [G3 (xprime)'];
end
% Quadratic term
if fitwhat(1) >= 2
  G2 = [G2 (xprime)'.^2];
  G3 = [G3 (xprime)'.^2];
end
% Cubic term
if fitwhat(1) == 3
  G3 = [G3 (xprime)'.^3];
end
% Periodic terms
if ~isempty(omega)
  % Angular frequency in radians/(rescaled day) of the periodic terms
  th_o= repmat(omega,nmonths,1)*2*pi.*repmat((xprime)',1,lomega);
  G1 = [G1 cos(th_o) sin(th_o)];
  G2 = [G2 cos(th_o) sin(th_o)];
  G3 = [G3 cos(th_o) sin(th_o)];
  % Create our specialterms G if we have it.  At the moment this is just
  % an additional periodic function, but in the future we could add something
  % else.
  if ~isnan(specialterms{1})
      % Strip off the previous periodic terms here and REPLACE with omegaspec
      omegaspec = [omega mu2/specialterms{3}];
      thspec= repmat(omegaspec,nmonths,1)*2*pi.*repmat((xprime)',1,length(omegaspec));
      Gspec1 = [G1(:,1:end-2*lomega) cos(thspec) sin(thspec)];
      Gspec2 = [G2(:,1:end-2*lomega) cos(thspec) sin(thspec)];
      Gspec3 = [G3(:,1:end-2*lomega) cos(thspec) sin(thspec)];
  end

end % end periodic
% Initilization Complete



%%%
% Solving
%%%

% Since each Slepian coefficient has different errors, each will have a 
% different weighting matrix.  Thus we loop over the coefficients.
for index=1:j 
    % If we have a priori error information, create a weighting matrix, and
    % change the G and d matrices to reflect this.  Since each coefficient
    % has its own weighting, we have to invert them separately.
    W = diag([1./givenerrors(:,index)]);
    d = slept(:,index);
    G1w = W*G1;
    G2w = W*G2;
    G3w = W*G3;
    dw = W*d;
 
    % This is in case you request a single special term to be looked at
    % in a special way
    if index == specialterms{1}
       Gspec1w = W*Gspec1;
       Gspec2w = W*Gspec2;
       Gspec3w = W*Gspec3;
       lomega = length(omegaspec);
       th=thspec;
       myomega=omegaspec;
    else
       lomega = length(omega);
       myomega=omega;
       th=th_o;
    end
    
    %%%
    % First order polynomial
    %%%
    
    % Do the fitting by minimizing least squares
    % First, the linear fit with periodics
    mL2_1 = (G1w'*G1w)\(G1w'*dw) ;
    % That was regular, but if there was a special one, substitute
    if index == specialterms{1}
       mL2_1 = (Gspec1w'*Gspec1w)\(Gspec1w'*dw) ;
    end

    % Use the model parameters to make the periodic amplitude in time
    startP = length(mL2_1) - 2*lomega + 1;
    amp1 = [mL2_1(startP:(startP+lomega-1)) mL2_1((startP+lomega):end)];
    amp1 = sqrt(amp1(:,1).^2 + amp1(:,2).^2);
    
    % Use the model parameters to make the periodic phase in time
    phase1 = [mL2_1(startP:(startP+lomega-1)) mL2_1((startP+lomega):end)];
    phase1 = atan2(phase1(:,1),phase1(:,2));

    % Assemble the estimated signal function, evaluated at 'thedates'
    % Start adding things
    fitfn1 = mL2_1(1) + mL2_1(2)*(xprime);

    % Add the sum over all the periodic components periodics
    fitfn1 = fitfn1 + ...
             sum(repmat(amp1,1,nmonths).*sin(th'+repmat(phase1,1,nmonths)),1);

    % For the special time when you want the periodic stuff taken out
    % fitfn1_nosine = mL2_1(1) + mL2_1(2)*(xprime);
    % resid1_nosine = d'-sum(repmat(amp1,1,nmonths).*sin(th'+repmat(phase1,1,nmonths)),1);
    
    % Compute the residual time series for this coefficient 
    resid1 = d - fitfn1';

    % Here's the definition of the residual at lm vs time
    ESTresid(:,index) = resid1;
    % Here's the definition of the signal at lm vs time
    ESTsignal(:,index) = fitfn1';
    
    % Do extra dates if you have them
    if moredates
       th_extras= repmat(myomega,length(extradatesprime),1)...
         *2*pi.*repmat((extradatesprime)',1,lomega);
       % Evaluate at the missing dates
       extravaluesfn1 = mL2_1(1) + mL2_1(2)*(extradatesprime) + ...
         sum(repmat(amp1,1,1).*sin(th_extras'+repmat(phase1,1,1)),1);
       extravalues(:,index) = extravaluesfn1';
    end
    
    % Get the residual sum of squares for later F tests
    rss1 = sum(resid1.^2);
    
    
    %%%
    % Second order polynomial
    %%%
    
    % Now repeat that all again with second order polynomial, if we want
    if fitwhat(1) >= 2  
      mL2_2 = (G2w'*G2w)\(G2w'*dw) ;
      if index == specialterms{1}
        mL2_2 = (Gspec2w'*Gspec2w)\(Gspec2w'*dw) ;
      end

      startP = length(mL2_2) - 2*lomega + 1;
      amp2 = [mL2_2(startP:(startP+lomega-1)) mL2_2((startP+lomega):end)];
      amp2 = sqrt(amp2(:,1).^2 + amp2(:,2).^2);
    
      phase2 = [mL2_2(startP:(startP+lomega-1)) mL2_2((startP+lomega):end)];
      phase2 = atan2(phase2(:,1),phase2(:,2));

      fitfn2 = mL2_2(1) + mL2_2(2)*(xprime) + mL2_2(3)*(xprime).^2;
      
      fitfn2 = fitfn2 + ...
             sum(repmat(amp2,1,nmonths).*sin(th'+repmat(phase2,1,nmonths)),1);

      % For the special time when you want the periodic stuff taken out
      % fitfn2_nosine = mL2_1(1) + mL2_1(2)*(xprime) + mL2_2(3)*(xprime).^2;
      % resid2_nosine = d'-sum(repmat(amp1,1,nmonths).*sin(th'+repmat(phase1,1,nmonths)),1);
    
      % Compute the residual time series for this coefficient 
      resid2 = d - fitfn2';
      
      % Do extra dates if you have them
      if moredates
         th_extras= repmat(myomega,length(extradatesprime),1)...
           *2*pi.*repmat((extradatesprime)',1,lomega);
         extravaluesfn2 = mL2_2(1) + mL2_2(2)*(extradatesprime) +...
           mL2_2(3)*(extradatesprime).^2 +...
           sum(repmat(amp1,1,1).*sin(th_extras'+repmat(phase1,1,1)),1);
      end
      
      % Get the residual sum of squares
      rss2 = sum(resid2.^2);
      % Calculate an F-score for this new fit
      fratioP2 = (rss1 - rss2)/1/(rss2/(length(slept(:,index))-length(mL2_2)));
      fscore = finv(.95,1,length(slept(:,index))-length(mL2_2));
      if fratioP2 > fscore
         P2ftest = 1;
         % We pass, so update the signal and residuals with this new fit
         ESTresid(:,index) = resid2;
         ESTsignal(:,index) = fitfn2';
         if moredates; extravalues(:,index) = extravaluesfn2'; end
      else
         P2ftest = 0;
      end
    
    end
    
    %%%
    % Third order polynomial
    %%%
    
    % Now repeat that all again with third order polynomial, if we want
    if fitwhat(1) >= 3
      mL2_3 = (G3w'*G3w)\(G3w'*dw) ;
      if index == specialterms{1}
        mL2_3 = (Gspec3w'*Gspec3w)\(Gspec3w'*dw) ;
      end

      startP = length(mL2_3) - 2*lomega + 1;
      amp3 = [mL2_3(startP:(startP+lomega-1)) mL2_3((startP+lomega):end)];
      amp3 = sqrt(amp3(:,1).^2 + amp3(:,2).^2);
    
      phase3 = [mL2_3(startP:(startP+lomega-1)) mL2_3((startP+lomega):end)];
      phase3 = atan2(phase3(:,1),phase3(:,2));

      fitfn3 = mL2_3(1) + mL2_3(2)*(xprime) ...
          + mL2_3(3)*(xprime).^2 + mL2_3(4)*(xprime).^3;
      
      fitfn3 = fitfn3 + ...
             sum(repmat(amp3,1,nmonths).*sin(th'+repmat(phase3,1,nmonths)),1);

      resid3 = d - fitfn3';

      % For the special time when you want the periodic stuff taken out
      % fitfn3_nosine = mL2_1(1) + mL2_1(2)*(xprime) +...
      %    mL2_2(3)*(xprime).^2 + mL2_3(4)*(xprime).^3;
      % resid3_nosine = d'-sum(repmat(amp1,1,nmonths).*sin(th'+repmat(phase1,1,nmonths)),1);
    
      % Do extra dates if you have them
      if moredates
         th_extras= repmat(myomega,length(extradatesprime),1)...
           *2*pi.*repmat((extradatesprime)',1,lomega);
         extravaluesfn3 = mL2_3(1) + mL2_3(2)*(extradatesprime) +...
           mL2_3(3)*(extradatesprime).^2 + mL2_3(4)*(extradatesprime).^3 +...
           sum(repmat(amp1,1,1).*sin(th_extras'+repmat(phase1,1,1)),1);
      end
      
      % Get the residual sum of squares
      rss3 = sum(resid3.^2);
      % Calculate an F-score for this new fit
      fratioP3 = (rss1 - rss3)/2/(rss3/(length(slept(:,index))-length(mL2_3)));
      fscore = finv(.95,1,length(slept(:,index))-length(mL2_3));
      if fratioP3 > fscore
         P3ftest = 1;
         % We pass, so update the signal and residuals with this new fit
         ESTresid(:,index) = resid3;
         ESTsignal(:,index) = fitfn3';
         if moredates; extravalues(:,index) = extravaluesfn3'; end
      else
         P3ftest = 0;
      end
      
    end
    
    % Some extra plotting for excessive verification
    if xver==1 && index <= 30
%       clf
%       subplot(211)
%       plot(thedates,d,'+')
%       hold on
%       plot(thedates,fitfn,'r');
%       errorbar(thedates,d,givenerrors(:,1),'b-')
% 
%       axis tight
%       datetick('x',28)
%       title(sprintf('components at alpha = %i',index))
%       legend('Data with given errors','Fitted signal')
%       
%       subplot(212)
%       DY=sqrt(mean(ESTresid(:,1).^2))*ones(size(fitfn))';
%       fill([thedates' ; flipud(thedates')],...
%            [ESTresid(:,1)+DY ; flipud(ESTresid(:,1)-DY)],...
%            [0.7 0.7 0.7]);
%       %errorfill(thedates,ESTresid(:,1),sqrt(mean(ESTresid(:,1).^2))*ones(size(fitfn)))
%       %errorbar(thedates,ESTresid(:,1),sqrt(mean(ESTresid(:,1).^2))*ones(size(fitfn)),'r-')
%       hold on
%       DY=givenerrors(:,1);
%       fill([thedates' ; flipud(thedates')],...
%            [ESTresid(:,1)+DY ; flipud(ESTresid(:,1)-DY)],...
%            [0.9 0.9 0.9]);
%       %errorbar(thedates,ESTresid(:,1),givenerrors(:,1),'b-')
%       title(sprintf('Residuals with error bands at alpha = %i',index))
%       legend('Estimated errors','Given errors')

%       axis tight
%       datetick('x',28)
        
        clf
        plot(thedates,d,'b-')
        hold on
        plot(thedates,ESTsignal(:,index),'r-')
        datetick('x',28)
        title(sprintf('alpha = %i',index))
        keyboard
    end
    
    % Make the matrix ftests
    ftests(index,:) = [0 P2ftest P3ftest];

end

% Collect output
varns={ESTsignal,ESTresid,ftests,extravalues};


%%%
% TOTAL COMBINED FITTING
%%%

% If we have the parameters for this localization, and we requested the
% total fit, then let's do that.

if nargout >= 5 && exist('CC') && exist('TH')
    
    % Get the residual covariance
    [Cab] = slepresid2cov(ESTresid);
    
    % Calculate the bandwdith for this basis
    L = sqrt(size(slept,2)) - 1;
    % This should be an integer
    if (floor(L)~=L);
        error('Something fishy about your L');
    end
    if iscell(TH)
       % Something like {'greenland' 0.5}
       XY=eval(sprintf('%s(%i,%f)',TH{1},10,TH{2}));
    else
       % Coordinates or a string, either works
       XY=TH;
    end
    % Calculate the Shannon number for this basis
    N=round((L+1)^2*spharea(XY));
    
    % Make the coefficients with reference to some mean
    % If they already are, then this won't matter
    sleptdelta = slept(1:nmonths,:) - repmat(mean(slept(1:nmonths,:),1),nmonths,1);

    
    % COMBINE

    % We want to take the Slepian functions and combine them to get total mass.
    % For signal, this means integrating the functions and adding them.  For
    % the error, this means using the error propogation equation, where we
    % compute (int)*(covar)*(int)'.  Since the slepcoffs are constants that 
    % just come forward, we can do the integration of the eigenfunctions 
    % first (and once for each function), then multiply by slepcoffs to 
    % get the monthly values.  This is much faster.

    [eigfunINT] = integratebasis(CC,TH,N);
    
    % Since Int should have units of (fn * m^2), need to go from fractional
    % sphere area to real area.  If the fn is surface density, this output is
    % in kilograms.  Then change the units from kg to Gt in METRIC tons
    eigfunINT = eigfunINT*4*pi*6370000^2/10^3/10^9;
    functionintegrals = eigfunINT;
    
    % Now multiply by the appropriate slepcoffs to get the months
    % This becomes alpha by months
    %functimeseries=repmat(eigfunINT',1,nmonths).*sleptdelta(:,1:N)';
    %functimeseries = sleptdelta(:,1:N)';

    % Here do the total sum of the data
    total=eigfunINT*sleptdelta(:,1:N)';

    % Get the error
    thevars = diag(Cab(1:N,1:N))';
    alphavar = eigfunINT.^2.*thevars;
    % Now the combined error with covariance
    alphavarall = eigfunINT*Cab(1:N,1:N)*eigfunINT';

    % FITTING

    % We have uniform estimated error, which will be different than the polyfit
    % estimated residuals because ours account for sinusoidal signals.  So
    % pass the new error to our function for replacement, so
    % that the fitting confidence intervals reflect that  
    
    [fit,delta,params,paramerrors] = timeseriesfit([thedates' total'],alphavarall,1,1);
    
    % Make a matrix for the line, and 95% confidence in the fit
    totalfit = [thedates' fit delta];

    % Make the parameters and error valid for a year
    totalparams = params(2)*365;
    totalparamerrors = paramerrors*365;

    % Collect the expanded output
    varns={ESTsignal,ESTresid,ftests,extravalues,...
       total,alphavarall,totalparams,totalparamerrors,totalfit,...
       functionintegrals,alphavar};    
end

varargout=varns(1:nargout);
