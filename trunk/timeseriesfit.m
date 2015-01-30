function varargout=timeseriesfit(timeseries,theerror,lines,polys)
% [thefits,thedeltas,params,paramerrors,ftests]=TIMESERIESFIT(timeseries,theerror,lines,polys)
%
% This function accepts a time series and does a choice of fitting, using
% and F-test to determine which of the chosen fits are significant.  All
% values are reported in the original x-coordinate units.
%
%
% INPUT:
%
% timeseries     A n-by-2 matrix where the first column are the X values
%                 and the second column are the Y values of the timeseries.
% theerror       The error on the timeseries data expressed as the 
%                 variance [default: determined from the residuals].
% lines          The number of linear fits you want to attempt (1-3) [default: 1].
%                 NOTE: currently this is disabled
% polys          The number of polynomials you want to try and fit [default: none].
%                 If you give 1, then nothing happens, as this is a linear fit.  
%                 2 is a quadratic fit and 3 is a cubic fit.
%            
% OUTPUT:
%
% thefits        An n-by-j matrix where n is the number of datapoints you
%                 had, and j is the number of fits you requested. i.e. each 
%                 column is a fit.  This begins with the line fits you
%                 requested and then progresses through the polynomials
%                 with increasing order.  The single order polynomial is
%                 not duplicated.
% thedeltas      A matrix similar to "thefits" which has the error envelope
%                 for the fits.  You plot this as Y+delta and Y-delta.
% params         The parameters of the fits (e.g. slope).  This is a 4-by-j
%                 matrix where j are the fits you wanted.  Zeros fill out 
%                 the unused parameters.  For example, if
%                 you want just a 2nd order polynomial, you will get 
%                 [intercept intercept; slope slope; 0 quadratic; 0 0]
% paramerrors    A matrix similar to params, that contains the two sigma
%                 errors on the parameters.  NOTE: currently we report 0
%                 for the error on the intercepts, but this could be added
%                 later.  This is because, for our needs, the intercept and
%                 its errors will be huge since the matlab dates are so far
%                 from the origin, and there is strong covariance between
%                 parameters in this case.
% ftests         An array, such as [0 1 1], on whether the fits you
%                 requested passed an F-test for significance.
%
%
%
% NOTES:  In this function we assume that we have uniform estimated error for 
%   each datapoint.  This is because polyfit does not accept error, let alone 
%   error that is different for each datapoint.  I believe another Matlab program
%   will accept error, but will not do a weighted least squares.  If you desire this
%   then go elsewhere.
%
% SEE ALSO:  SLEPT2RESID, LOCALIZATION2FIT
%
% Last modified by charig-at-princeton.edu on 1/5/2012

defval('lines',1)
defval('polys',3)
defval('P3ftest',[]);


%%%
% INITIALIZE
%%%

%load('testdata.mat')
%defval('theerror',0);

% We may have uniform estimated error, which will be different than the polyfit
% estimated residuals because ours account for sinusoidal signals.  So use
% polyfit here, but if we have it, replace the residual norm with our error estimate, so
% that the fitting confidence intervals reflect that.


% Split the timeseries for ease
theX = timeseries(:,1);
theY = timeseries(:,2);
%theX = thedates(:);
%theY = total(:);
%theerror = alphavarall;
% For reference, we will assume the X data is time in Matlab date format (decimal 
% days from some origin) and the Y data is in mass.

%%%
% FITTING
%%%

% The first thing we do is a single linear fit.

[PhatL1,SL1,MUL1] = polyfit(theX,theY,1);
% These are the rescaled values (slope in per day)
PL1(1) = PhatL1(1)/MUL1(2);
PL1(2) = PhatL1(2) - PL1(1)*MUL1(1);
% Replace it here, if you have your own error
onelinenormr = SL1.normr;
if ~isempty(theerror)
   SL1.normr = sqrt(SL1.df*theerror);
end
% Get the fit confidence intervals
[L1,L1delta] = polyconf(PhatL1,theX,SL1,'mu',MUL1,'predopt','curve');
% Make the covar matrix in order to calculate the slope error
% Now use S.normr and S.df since we made a change above
fitcovarL1 = (inv(SL1.R)*inv(SL1.R)')*SL1.normr^2/SL1.df;
% Back to original coordinates for the slope
slopeL1var = fitcovarL1(1)/MUL1(2)^2;
% As two sigma error in mass/day
slopeL1twosigma = sqrt(slopeL1var)*tinv(.975,SL1.df);
% The slope in mass/year
%singleslope = PL1(1)*365;

% Gather
thefits = L1;
thedeltas = L1delta;
params = [PL1(2); PL1(1); 0; 0];
paramerrors = [0; slopeL1twosigma; 0; 0];

%%%
% Additional Lines
%%%

% % Now maybe do more than one linear fit for this series.
% myftest = zeros(108,1);
% myslopesok = zeros(108,1);
% 
% if lines == 2 | lines == 3
%    % Do 2 line fits by finding the best place to split the timeseries
%    for i = 3:(length(theX)-3)
%       % Do a linear fit, and save the best one in an normr sense
%       % Then later we will use the df to find if this new fit is significant
%       % And also check whether the two lines have the same slope within error
%       [Phat1,S1,MU1] = polyfit(theX(1:i),theY(1:i),1);
%       [Phat2,S2,MU2] = polyfit(theX(i+1:end),theY(i+1:end),1);
% 
%       % What is the new combined normr
%       newnormr = sqrt(S1.normr^2 + S2.normr^2);
%       
%       % These are the rescaled values (slope in per day)
%       Ptot1(1) = Phat1(1)/MU1(2);
%       Ptot1(2) = Phat1(2) - Ptot1(1)*MU1(1);
%       Ptot2(1) = Phat2(1)/MU2(2);
%       Ptot2(2) = Phat2(2) - Ptot2(1)*MU2(1);
%       % Replace it here
%       S1.normr = sqrt(S1.df*theerror);
%       S2.normr = sqrt(S2.df*theerror);
%       % Get the fit confidence intervals
%       [Y1,delta1] = polyconf(Phat1,theX(1:i),S1,'mu',MU1,'predopt','curve');
%       [Y2,delta2] = polyconf(Phat2,theX(i+1:end),S2,'mu',MU2,'predopt','curve');
%       % Make the covar matrix in order to calculate the slope error
%       % Now use S.normr and S.df since we made a change above
%       fitcovar1 = (inv(S1.R)*inv(S1.R)')*S1.normr^2/S1.df;
%       fitcovar2 = (inv(S2.R)*inv(S2.R)')*S2.normr^2/S2.df;
%       % Back to original coordinates
%       slopevar1 = fitcovar1(1)/MU1(2)^2;
%       slopevar2 = fitcovar2(1)/MU2(2)^2;
%       % As two sigma error in mass/year
%       slopetwosigma1 = sqrt(slopevar1)*365*1.99;
%       slopetwosigma2 = sqrt(slopevar2)*365*1.99;
%       % The slopes
%       twoslope1 = Ptot1(1)*365;
%       twoslope2 = Ptot2(1)*365;
%       % Make a matrix for the line
%       %myfit = [thedates' Y'];
%       % Make a matrix for the 95% error lines
%       %myfiterror = [delta'];
% 
%       % Do the f test to see if this break is significant
%       fscore = (onelinenormr^2 - newnormr^2)/3/(newnormr^2/103)
%       if fscore > 3.254
%          myftest(i) = 1;
%       end
%       if (twoslope1 - twoslope2) > (slopetwosigma1+slopetwosigma2)
%          myslopesok(i) = 1;
%       end
%       if i==12
%       twoslope1
% slopetwosigma1
% twoslope2
% slopetwosigma2
% newnormr
% newnormr^2
%       end
% 
% 
% 
%    end
% 
% end


%%%
% Additional Polynomials
%%%

if polys == 2 || polys == 3
    % fit a quadratic polynomial
    % Here we do a 2nd order polynomial fit
    [PhatP2,SP2,MUP2] = polyfit(theX,theY,2);
    % These are the rescaled values (in per day)
    PP2(1) = PhatP2(1)/MUP2(2)^2;
    PP2(2) = PhatP2(2)/MUP2(2);
    PP2(3) = PhatP2(3) - PP2(2)*MUP2(1);
    % Get the fit confidence intervals
    [P2,P2delta] = polyconf(PhatP2,theX,SP2,'mu',MUP2,'predopt','curve');
    % Make the covar matrix in order to calculate the slope error
    % Now use S.normr and S.df since we made a change above
    fitcovarP2 = (inv(SP2.R)*inv(SP2.R)')*SP2.normr^2/SP2.df;
    % Back to original coordinates
    accelvarP2 = fitcovarP2(1,1)/MUP2(2)^4;
    slopevarP2 = fitcovarP2(2,2)/MUP2(2)^2;
    % As two sigma error in mass/year
    accelP2twosigma = sqrt(accelvarP2)*tinv(.975,SP2.df);
    slopeP2twosigma = sqrt(slopevarP2)*tinv(.975,SP2.df);
    % The slope and accel
    accelP2 = PP2(1)*365*365;
    slopeP2 = PP2(2)*365;

    % Calculate an F-score for this new fit
    fratioP2 = (onelinenormr^2 - SP2.normr^2)/1/(SP2.normr^2/SP2.df);
    fscore = finv(.95,1,SP2.df);
    if fratioP2 > fscore
       P2ftest = 1;
    else
       P2ftest = 0;
    end

    % Gather
    thefits = [thefits P2];
    thedeltas = [thedeltas P2delta];
    params = [params [PP2(3); PP2(2); PP2(1); 0]];
    paramerrors = [paramerrors [0; slopeP2twosigma; accelP2twosigma; 0]];

    
    if polys == 3
        % fit a cubic polynomial
        % Here we do a 3nd order polynomial fit
        [PhatP3,SP3,MUP3] = polyfit(theX,theY,3);
        % These are the rescaled values (in per day)
        PP3(1) = PhatP3(1)/MUP3(2)^3;
        PP3(2) = PhatP3(2)/MUP3(2)^2;
        PP3(3) = PhatP3(3)/MUP3(2);
        PP3(4) = PhatP3(3) - PP3(3)*MUP3(1);
        % Replace it here
        %onelinenormr = S.normr;
        %S.normr = sqrt(S.df*theerror);
        % Get the fit confidence intervals
        [P3,P3delta] = polyconf(PhatP3,theX,SP3,'mu',MUP3,'predopt','curve');
        % Make the covar matrix in order to calculate the slope error
        % Now use S.normr and S.df since we made a change above
        fitcovarP3 = (inv(SP3.R)*inv(SP3.R)')*SP3.normr^2/SP3.df;
        % Back to original coordinates
        cubicvarP3 = fitcovarP3(1,1)/MUP3(2)^8;
        accelvarP3 = fitcovarP3(2,2)/MUP3(2)^4;
        slopevarP3 = fitcovarP3(3,3)/MUP3(2)^2;
        % As two sigma error in mass/year
        cubicP3twosigma = sqrt(cubicvarP3)*tinv(.975,SP2.df);
        accelP3twosigma = sqrt(accelvarP3)*tinv(.975,SP2.df);
        slopeP3twosigma = sqrt(slopevarP3)*tinv(.975,SP2.df);
        % The cubic, accel, and slope
        cubicP3 = PP3(1)*365*365*365;
        accelP3 = PP3(2)*365*365;
        slopeP3 = PP2(2)*365;
        
        % Calculate an F-score for this new fit
        fratioP3 = (SP2.normr^2 - SP3.normr^2)/1/(SP3.normr^2/SP3.df);
        fscore = finv(.95,1,SP3.df);
        if fratioP3 > fscore
           P3ftest = 1;
        else
           P3ftest = 0;
        end
        
        % Gather
        thefits = [thefits P3];
        thedeltas = [thedeltas P3delta];
        params = [params [PP3(4); PP3(3); PP3(2); PP3(1)]];
        paramerrors = [paramerrors [0; slopeP3twosigma; accelP3twosigma; cubicP3twosigma]];

        
    end % end polys==3

end % end polys==2



%%%
% OUTPUT
%%%
if polys == 2 || polys == 3
   ftests = [0 P2ftest P3ftest];
   varns={thefits,thedeltas,params,paramerrors,ftests};
else
   varns={thefits,thedeltas,params,paramerrors};
end
varargout=varns(1:nargout);



