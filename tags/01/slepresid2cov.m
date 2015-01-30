function varargout=slepresid2cov(ESTresid,J,imonths)
% [Cab,Cabr,Cabd]=SLEPRESID2COV(ESTresid,J,imonths)
%
% Computes a Slepian covariance matrix from a matrix with time
% series of residuals determined for each of the Slepian coefficients.
%
% INPUT:
%
% ESTresid    Residual time series for each pair of cos/sin coefficients
%             as determined from PLMT2RESID [default is to compute it]
% J           Bandwidth of the resulting covariance matrix [default: all]
% imonths     Index vector (such as [1:10] or [12:24]) of months to be 
%              considered for calculation [default: all]
%              noting that all of the months were likely part of the 
%              fitting procedure
%
% OUTPUT:
%
% Cab      The spectral covariance matrix
% Cabr     The scaled spectral covariance matrix (the correlation)
% Cabd     The diagonal of the covariance matrix
%
% EXAMPLE:
%
% slepresid2cov('demo1') % Also makes a plot
%
% Last modified by charig-at-princeton.edu, 06/27/2012

defval('ESTresid','slept2resid')
if isstr(ESTresid) && ~strcmp(ESTresid(1:4),'demo')
  % Evaluate the specified expression
  [ESTsignal,ESTresid] = eval(slept);
end

if ~isstr(ESTresid)

  %defval('L',ESTresid(1,end,1))
  defval('imonths',1:size(ESTresid,1))
  defval('xver',1)

  if max(imonths)~=size(ESTresid,1)
    warning(...
        'You are using a subset of residuals determined from the full set')
  end

  % Use your imonths clip (imonths should be a vector as in [1:20] or [12:24])
  ESTresid = ESTresid(imonths,:);
  
  % Spectral covariance of this "noise" as in, unfittable by SLEPT2RESID
  Cab=cov(ESTresid);

  if nargout>=2
    % Scaled version, note this is CORRCOEFF exactly
    Cabd=diag(sqrt(Cab));
    Cabr=Cab./[Cabd*Cabd'];
    Cabr(isnan(Cabr))=0;
  else
    Cabr=NaN;
    Cabd=NaN;
  end

  % Output
  varns={Cab,Cabr,Cabd};
  varargout=varns(1:nargout);

elseif strcmp(ESTresid,'demo1')
  % Could use the default matrix from SLEPRESID2COV, but better to specify
  % our parameters for the basis

  % Get GRACE data in a basis
  L=30;
  TH='amazon';
  [slepcoffs,calerrors,thedates,G,CC,V] = grace2slept('CSR',TH,0,L);
  % Remove a fit
  [ESTresid,thedates,ESTsignal,rmset,rmsst,varet,varst]=...
      slept2resid(slepcoffs,thedates,[1 1 365.0 181.0],calerrors);
  % Covariance matrix
  [Cab,Cabr,Cabd]=slepresid2cov(ESTresid);
  % Find the Shannon number
  N=length(Cab)*spharea(TH);
  
  clf
  % Plot the matrix, using functions up to 2*N
  Cabr=Cabr(1:round(2*N),1:round(2*N));
  crange=halverange(Cabr,75);
  imagefnan([0 0],[1 1],setnans(Cabr,100),[],crange)
  axis ij
  degshow=sort([1 round(N) round(2*N)]);
  tick=degshow/length(Cabr);
  degshow={'1' 'N' '2N'};
  set(gca,'XTick',tick, 'XTickLabel',degshow,'YTick',tick,'YTickLabel',degshow);
  xl(1) = xlabel('spherical harmonic degree l''');
  yl(1) = ylabel('spherical harmonic degree l');
  longticks(gca,2)
  [cb,xcb] = addcb('vert',crange,crange,'kelicol',0.5);
  set(cb,'yaxisl','r')
  watis='noise correlation matrix';
  dateform='mmmm yyyy';
  axes(cb)
  set(xcb,'string',sprintf('%s from %i months\n between %s and %s',...
                           watis,length(thedates),...
                           datestr(datevec(thedates(1)),dateform),...
                           datestr(datevec(thedates(end)),dateform)));
  shrink(cb,1.25,1)
  figdisp
end

