function varargout=resid2plot(ESTresid,thedates,ESTsignal,l,m,varet,calerrors)
% [ESTresid,thedates,ESTsignal,rmset,rmsst,varet,varst]
%             =RESID2PLOT(plmt,thedates,fitwhat,givenerrors)
%
% INPUT: 
% ESTresid    Residual time series for each pair of cos/sin coefficients
%             [nmonths x Ldata x 4] with a hyper-lmcosi format
%             ESTresid(:,:,1) the degree
%             ESTresid(:,:,2) the order
%             ESTresid(:,:,3) the cosine coefficient
%             ESTresid(:,:,4) the sine coefficient
% thedates    The dates that belong to these time series 
%               (in Matlab DATENUM format)
% l           The degree you want plotted [default: random]
% m           The order you want plotted  [default: random]
% ESTsignal   The estimated signal function for each cos/sin
%             coefficient evaluated at those months in the same format.  
%             This is optional if you want it plotted as well. 
% varet       The variance of the estimated residual coefficients.  If you
%               have these, they are plotted.  If not, they are estimated 
%               from the time series.
% calerrors   The given calibrated errors of the spherical harmonic 
%               time series, if you have them.
%

defval('ESTresid',[])

if ~isstr(ESTresid)

% Set default values
if nargin < 2
    % If you don't have 2 arguments, use a default for the first 3.
    % Otherwise, if you give it exactly 2 arguments, we assume you don't
    % want ESTsignal plotted, and it is not used.
    [ESTresid,thedates,ESTsignal]=plmt2resid;
end

% If no l or m, choose a random one from ESTresid, but not degree 0 or 1
therow=randi(length(ESTresid(1,:,1))-4)+4;
defval('l',ESTresid(1,therow,1));
defval('m',ESTresid(1,therow,2));

[~,~,~,~,~,~,~,~,~,~,demin]=addmon(l,m);
therow=demin(end);
% Make default errors
defval('varet',[squeeze(ESTresid(1,:,1:2))...
                squeeze(var(ESTresid(:,:,3:4),1)) ]);
defval('xver',0)  

% Take only what you want
residplot = squeeze(ESTresid(:,therow,:));
imonths = length(residplot);
varet = ones(imonths,2).*repmat(sqrt(varet(therow,3:4)),imonths,1);

% Round to 1 January of the first available and next last available year 
bounds=[datevec(thedates(1)); datevec(thedates(end))];
x0 = datenum(bounds(1)  ,1,1);
x1 = datenum(bounds(2)+1,1,1)-1;

% Now make a plot
subplot(211)
errorbar(thedates,residplot(:,3),varet(:,1),'r-')
hold on
mylegendtext = {'Estimated Std. Dev.'};
if exist('calerrors')
    % Plot the given errors (assumed to be smaller) over the estimated errors.
    errorbar(thedates,residplot(:,3),calerrors(:,therow,3),'b-')
    mylegendtext = [mylegendtext 'GRACE Calibrated Errors'];
end
if exist('ESTsignal')
    % Quickly remove the mean of signal so we can plot residuals and signal
    % together in one panel
    signalplot = squeeze(ESTsignal(:,therow,3:4))...
                 -repmat( squeeze(mean(ESTsignal(:,therow,3:4),1))',imonths,1);
    % Plot the original data (residual + signal)
    plot(thedates,residplot(:,3)+signalplot(:,1),'b-');
    mylegendtext = [mylegendtext 'Original data'];

end
plot([x0 x1],[0 0],'k-')
legend(mylegendtext)
datetick('x',28)
xl(1) = xlabel('Time');
yl(1) = ylabel('Geoid anomaly (m)');
longticks(gca,2)
axis tight
title(sprintf('cosine component at l = %i and m = %i',l,m))

subplot(212)
errorbar(thedates,residplot(:,4),varet(:,1),'r-')
hold on
mylegendtext = {'Estimated Std. Dev.'};

if exist('calerrors')
    % Plot the given errors (assumed to be smaller) over the estimated errors.
    errorbar(thedates,residplot(:,4),calerrors(:,therow,4),'b-')
    mylegendtext = [mylegendtext; 'GRACE Calibrated Errors'];
end
if exist('ESTsignal')
    % Plot the original data (residual + signal)
    plot(thedates,residplot(:,4)+signalplot(:,2),'b-');
    mylegendtext = [mylegendtext; 'Original data'];
end

plot([x0 x1],[0 0],'k-')
legend(mylegendtext)
datetick('x',28)
xl(1) = xlabel('Time');
yl(1) = ylabel('Geoid anomaly (m)');
longticks(gca,2)
axis tight
title(sprintf('sine component at l = %i and m = %i',l,m))



% Collect output
varns={1};
varargout=varns(1:nargout);

elseif strcmp(ESTresid,'demo1')
    % Get the original data
    [potcoffs,calerrors,thedates]=grace2plmt('CSR');
    % Get the fitted results
    [ESTresid,thedates,ESTsignal,~,~,varet]=plmt2resid(potcoffs(:,:,1:4),thedates,[1 1 181.0 365.0]);
    % Use both to make a plot
    clf
    resid2plot(ESTresid,thedates,ESTsignal,20,20,varet,calerrors);
end
