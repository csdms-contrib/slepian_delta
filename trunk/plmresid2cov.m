function varargout=plmresid2cov(ESTresid,L,imonths)
% [Clmlmlp,Clmlmpr,Clmlmpd,EL,EM]=PLMRESID2COV(ESTresid,L,imonths)
%
% Computes a spherical-harmonic covariance matrix from a matrix with time
% series of residuals determined for each of the harmonic coefficients.
%
% INPUT:
%
% ESTresid    Residual time series for each pair of cos/sin coefficients
%             as determined from PLMT2RESID [default is to compute it]
% L           Bandwidth of the resulting covariance matrix [default: all]
% imonths     Index of months to be considered for calculation [default: all]
%             noting that all of the months were part of the fitting procedure
%
% OUTPUT:
%
% Clmlmp      The spectral covariance matrix
% Clmlmpr     The scaled spectral covariance matrix (the correlation)
% Clmlmpd     The diagonal of the covariance matrix
% EL, EM      The spherical harmonic degree and order listing
%
% EXAMPLE:
%
% resid2cov('demo1') % Also makes a plot
%
% Last modified by fjsimons-at-alum.mit.edu, 06/26/2011
% Last modified by charig-at-princeton.edu, 06/26/2011

defval('ESTresid',plmt2resid)
warning('off','MATLAB:divideByZero');

if ~isstr(ESTresid)

  defval('L',ESTresid(1,end,1))
  defval('imonths',1:size(ESTresid,1))
  defval('xver',1)

  if max(imonths)~=size(ESTresid,1)
    warning(...
        'You are using a subset of residuals determined from the full set')
  end

  % Unwrap the residual coefficients so they can be run through the cov() function
  [~,~,~,~,~,~,~,~,~,ronm]=addmon(L);
  % The target indexing is like the below
  [EM,EL]=addmout(L);
  Lup=addmup(L);
  % Test that this is indeed the case for a random month
  if xver==1
    randm=ceil(rand*size(ESTresid,1));
    difer(EL-indeks(repmat(ESTresid(1,1:Lup,1)',1,2),ronm),[],[],NaN);
    difer(abs(EM)-indeks(repmat(ESTresid(1,1:Lup,2)',1,2),ronm),[],[],NaN);
  end

  % Initialize
  ESTresid2=nan(length(imonths),addmoff(L));
  % Reorder
  for index=1:length(imonths)
    ESTresid2(index,:)=indeks(ESTresid(imonths(index),1:Lup,3:4),ronm)';
  end

  % Spectral covariance of this "noise" as in, unfittable by PLMT2RESID
  Clmlmp=cov(ESTresid2);

  if nargout>=2
    % Scaled version, note this is CORRCOEFF exactly
    Clmlmpd=diag(sqrt(Clmlmp));
    Clmlmpr=Clmlmp./[Clmlmpd*Clmlmpd'];
    Clmlmpr(isnan(Clmlmpr))=0;
  else
    Clmlmpr=NaN;
    Clmlmpd=NaN;
  end

  % Output
  varns={Clmlmp,Clmlmpr,Clmlmpd,EL,EM};
  varargout=varns(1:nargout);

elseif strcmp(ESTresid,'demo1')
  % Get the default matrix
  [Clmlmp,Clmlmpr,Clmlmpd,EL,EM]=resid2cov([],20);
  % Overkill in a way, but get the dates also, make it fast
  [~,thedates]=plmt2resid([]);
  % Maximum degree of the expansion
  L=addmoff(size(Clmlmpr,1),'r');
  % Block sort the correlation matrix?
  [EM2,EL2,mz,blkm,dblk]=addmout(L);
  % Check this is what it is
  difer(EM-EM2); difer(EL-EL2)
  % Block sort if you want
  % Clmlmpr=Clmlmpr(blkm,blkm);
  clf
  % Plot the matrix
  crange=halverange(Clmlmpr,75);
  imagefnan([0 0],[1 1],setnans(Clmlmpr,100),[],crange)
  
  axis ij
  degshow=sort([2 L/2 3/4*L L]);
  tick=(addmoff(degshow-1)+1)/(L+1)^2;
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

