function varargout=Clmlmp2Crrp(Clmlmp,rnot,degres,VARoption)
% [Covrr,Covrr_scaled,Covrr_shaped,rnot,expco]=CLMLMP2CRRP(Clmlmp,rnot,degres)
% [r,lon,lat]=CLMLMP2CRRP(Clmlmp,[],degres,VARoption)
%
% Takes a spherical harmonic covariance matrix and calculates the
% spatial covariance. (map = Yrplm*Clmlmp*Yrlm')  
% OR in the case of an empty 'rnot', makes a field of spatial variance
%
% INPUT:
%
% Clmlmp    Spherical harmonic covariance matrix, assumed to be ordered as
%             in GLMALPHA (not KERNELC)
% rnot      The point(s) with which we want to find spatial covariance.
%             - If rnot is empty(i.e. []) then do whole sphere variance [default]
%             - Otherwise, rnot is a n-by-2 matrix of 
%               [lon colat; lon2 colat2; etc]
% degres    The degree resolution you want. [default: Nyquist]
% option    The method used to calculate the global variance
%             - 1: Calculate all of Crrp and select diagonal (memory
%                  inefficient)
%             - 2: Directly calculate just the diagonal of Crrp [default]
%             - 3: Loop over the spatial points and calculate variance at
%                  each point (slow)
%
% OUTPUT:
%
% Covrr           The spatial covariance, where dimensions 
%                   are [#locations x WholeSphereGrid]
% Covrr_scaled    The spatial covariance, where each row is scaled by the
%                   variance of the rnot point
% Covrr_shaped    The scaled spatial covariance, reshaped into size of the
%                   phi,theta meshgrid, as a cell array, for easy plotting 
%                   through imagefnan
% rnot            Your points back to you
% lmcosi          What is actually expanded to be plotted
%
%          -- OR --
%
% r,lon,lat       When you want just variance on the sphere, these are the
%                   same output format as from PLM2XYZ
%
% EXAMPLE:
%
% Clmlmp2Crrp('demo1') Check the 3 variance methods against eachother 
% Clmlmp2Crrp('demo2') Make a figure of the global variance up to L=60
% Clmlmp2Crrp('demo3') Same as demo2 but centered over Greenland
% Clmlmp2Crrp('demo4') Make a figure of the spatial covariance of a point
%                      centered in southern Greenland
% Clmlmp2Crrp('demo5') As demo4 but neglecting spectral off-diagonals
%
% SEE ALSO:
%
% NOTES: Suggest plotting Covrr_shaped with 
%          imagef([0 90],[360 -90],Covrr_shaped{1});
%
% Last modified by charig-at-princeton.edu 10/24/2011

defval('Clmlmp',plmresid2cov)

if ~isstr(Clmlmp)
  
  % Determine the bandwidth of our data.  Assume Clmlmp starts at l=0.
  [n,m]=size(Clmlmp);
  L=addmoff(n,'r');

  % Default resolution is the Nyquist degree; return equal sampling in
  % longitude and latitude; sqrt(L*(L+1)) is equivalent wavelength
  degN=180/sqrt(L*(L+1));
  defval('degres',degN);
  defval('VARoption',2);
  defval('xver',0);
  defval('rnot',[]);
  
  % Make a grid
  % Default grid is all of the planet
  defval('c11cmn',[0 90 360 -90]);
  % Build in maxima to save memory
  defval('latmax',Inf); 
  defval('lmax',999);
  if degres>degN
    disp('Clmlmp2Crrp: You can do better! Ask for more spatial resolution')
    disp(sprintf('Spatial sampling ALLOWED: %8.3f ; REQUESTED: %6.3f',...
                  degN,degres))
  end
  % The number of longitude and latitude grid points that will be computed
  nlon=min(ceil([c11cmn(3)-c11cmn(1)]/degres+1),latmax);
  nlat=min(ceil([c11cmn(2)-c11cmn(4)]/degres+1),2*latmax+1);
    
  % Longitude grid vector in radians
  phi=linspace(c11cmn(1)*pi/180,c11cmn(3)*pi/180,nlon);
  % Colatitude grid vector in radians
  theta=linspace([90-c11cmn(2)]*pi/180,[90-c11cmn(4)]*pi/180,nlat);
    
  % Make these into longer vectors for each point, in degrees
  phi = gamini(phi*180/pi,nlat);
  theta = repmat(theta*180/pi,1,nlon);

  % Handle the locations
  [c,d]=size(rnot);
  locations=zeros(1,c);
  if c~= 0 % We have some points where we want covariance
      
    for i=1:c
      foundit = [find(single(phi)==rnot(i,1) & single(theta)==rnot(i,2))];
      if isempty(foundit)
        error('Requested rnot does not work with degres')
      else
        locations(i) = foundit;
      end
    end        

    % The formula we wish to calculate is:
    % Yrprimelm*Clmlmp*Yrlm
    % To be quick, we can calculate just the first multiplication with the
    % rprimes of interest (rnot), and then expand these coefficients by
    % using plotplm.  This is faster because plotplm uses save files for
    % the Legendre expansions while ylm does not.
    
    % Get some Yrlms, just for our rnot points
    [Yrlm,~,~,dems,dels] = ylm([0 L],[],theta(locations)*pi/180,phi(locations)*pi/180,[],[],[],1);
    % Adjust Yrlm for phase and normalization, Make sure we're on the same page
    [EM,EL,~,]=addmout(L); difer(dems-EM,[],[],NaN); difer(dels-EL,[],[],NaN)
    Yrlm=[Yrlm.*repmat((-1).^EM,1,size(Yrlm,2))]*sqrt(4*pi);
  
    % Calculate our covariance for these rnot points
    expco = Yrlm'*Clmlmp;
    [~,~,~,lmcosi,~,~,~,~,~,ronm]=addmon(L);

    if xver==1
      % You can verify we did this correctly by doing the whole calculation
      % Get some Yrlms, just for our rnot points
      [Yrlm2,~,~,dems,dels] = ylm([0 L],[],theta*pi/180,phi*pi/180,[],[],[],1);
      % Now check the phase and the normalizations
      Yrlm2=[Yrlm2.*repmat((-1).^EM,1,size(Yrlm2,2))]*sqrt(4*pi);
      % Make a Y matrix just for the points of interest
      Yrprimelm2 = Yrlm2(:,locations)';
      % Do eet
      Covrr2 = Yrprimelm2*Clmlmp*Yrlm2;
      Covrr2=reshape(Covrr2(1,:),nlat,nlon); 
      % Now we need to look at it.
    end

    % Expand, Rescale, and reshape
    for want=1:length(locations)
      % Perform the expansion of [Yrprimelm*Clmlmp]
      lmcosi(ronm+size(lmcosi,1)*2)=expco(want,:);
      Covrr_SH{want} = lmcosi;
      [Covrrx,~,~,lon2,lat2]=plotplm(lmcosi,[],[],4,1);
      close;
      % Now you can use PLM2ROT to rotate the lmcosi if you choose
      Covrr(want,:) = reshape(Covrrx,1,nlat*nlon);

      Covrr_scaled(want,:) = Covrr(want,:)/Covrr(want,locations(want));
      Covrr_shaped{want} = reshape(Covrr_scaled(want,:),nlat,nlon);
    end
    
    % Collect output
    varns={Covrr,Covrr_scaled,Covrr_shaped,rnot,Covrr_SH};
    varargout=varns(1:nargout);

  else % We just want global variance
      
    % Get some Yrlms
    [Yrlm,~,~,dems,dels] = ylm([0 L],[],theta*pi/180,phi*pi/180,[],[],[],1);
    % Adjust Yrlm for phase and normalization, Make sure we're on the same page
    [EM,EL,~,]=addmout(L); difer(dems-EM,[],[],NaN); difer(dels-EL,[],[],NaN)
    Yrlm=[Yrlm.*repmat((-1).^EM,1,size(Yrlm,2))]*sqrt(4*pi);
  
    % We calculate the global variance taking into account the spectral
    % covariance terms.  Note this is not the diagonal elements of the spectral
    % covariance matrix, but instead should be the diagonal elements of the
    % spatial covariance matrix.  We have 3 ways to calculate this:
    
    if VARoption == 1
      % Option 1: Directly calculate the whole matrix as above, and use the
      % diagonal elements.  This is memory inefficient, but computationally
      % fast.
      spacepoints = ((180/degres+1)*(360/degres+1));
      if (spacepoints*addmoff(L,'a')) > 2e7 % degres < 10 when L=60
        error(['Your degree resolution will likely run afoul of the' ...
               ' memory limitations.  Consider lower resolution or a different option']);
      end
      % Make the second Y matrix for the whole globe
      Yrprimelm = Yrlm';
      % Do the calculation all at once
      Covrr = Yrprimelm*Clmlmp*Yrlm;
      r = diag(Covrr);

    elseif VARoption == 2
      % Option 2: Calculate just the diagonal elements of the spatial
      % covariance matrix directly.  
      
      % We are calculating the square of the expansion, which consists 
      % of products of Ylm and Ylmp, and elements of the covariance 
      % matrix.  We do this by performing half the calculation, and 
      % then use a dot-multiply to compute what are the diagonal 
      % elements.  You can verify this if you like.
      part1 = Clmlmp*Yrlm;
      r = sum(Yrlm.*part1,1)';
      
      if xver ==1
        Yrprimelm = Yrlm';
        Covrr = Yrprimelm*Clmlmp*Yrlm;
        difer(r-diag(Covrr))
      end
      
      
    elseif VARoption == 3
      % Option 3: Loop over each spatial point and calculate
      % Y(r)*Clmlmp*Y(r') for each point in one go.  This is
      % computationally inefficient, but uses much less memory.
      h=waitbar(0,'Clmlmp2Crrp: Loop over spatial points to get variance');
      Yrprimelm = Yrlm';
      lenphi = length(phi);
      for f = 1:lenphi
        Covrr(f) = Yrprimelm(f,:)*Clmlmp*Yrlm(:,f);
        waitbar((f/lenphi),h);
      end
      r = Covrr';
      close(h)

    end % end VARoption
      
    % Respahe to give matrix output
    r     = reshape(r,length(unique(theta)),length(unique(phi)));
    phi   = reshape(phi,length(unique(theta)),length(unique(phi)));
    theta = reshape(theta,length(unique(theta)),length(unique(phi)));
    
    % Collect output
    varns={r,phi,theta};
    varargout=varns(1:nargout);
      
      
  end % end covariance at points or global variance?

    
elseif strcmp(Clmlmp,'demo1')
  % Test the efficiency and accuracy of the three methods
  L=30;
  [Clmlmp]=plmresid2cov([],L);
  
  disp('Method 1:')
  tic
  [r1,lon1,lat1]=Clmlmp2Crrp(Clmlmp,[],[],1);
  toc
  disp('Method 2:')
  tic
  [r2,lon2,lat2]=Clmlmp2Crrp(Clmlmp,[],[],2);
  toc
  disp('Method 3:')
  tic
  [r3,lon3,lat3]=Clmlmp2Crrp(Clmlmp,[],[],3);
  toc
  disp('Check 1 vs. 2')
  difer(r1-r2)
  disp('Check 1 vs. 3')
  difer(r1-r3)
  
  
elseif strcmp(Clmlmp,'demo2')
  % Lets make a figure of the global variance
  L=60;
  Clmlmp=plmresid2cov([],L);
  [r,lon,lat]=Clmlmp2Crrp(Clmlmp,[],[]);
  figure
  imagef([0 90],[360 -90],r);
  % caxis(ah1(want),[-2.5 2.5]);
  plotcont
  axis([0 360 -90 90])
  colorbar
  longticks(gca,2)
  t=title('Spatial Variance of noise');
  % saveas(gcf,['SpatialVariance_L' num2str(L)],'epsc');
  
  % Also make this figure in mm water equivalent
  % Mean Earth radius
  a=fralmanac('a_EGM96','Earth');
  % Change variance to standard deviations in meters.
  r = a*sqrt(r);
  % Convert r to spherical harmonics
  [lmcosi] = xyz2plm(r,60);
  lmcosiSD=plm2pot_harig(lmcosi,[],[],[],4);
  [r2,lon2,lat2,Plm2]=plm2xyz(lmcosiSD);  
  figure
  imagef([0 90],[360 -90],r2);
  % caxis(ah1(want),[-2.5 2.5]);
  plotcont
  axis([0 360 -90 90])
  colorbar
  longticks(gca,2)
  t=title('Standard deviation of noise (mm water equivalent)');
  
elseif strcmp(Clmlmp,'demo3')
  % Lets make a figure of the error centered over Greenland, expressed 
  % in mm of water equivalent
  L=30;
  Clmlmp=plmresid2cov([],L);
  [r,lon,lat]=Clmlmp2Crrp(Clmlmp,[],1);
  figure
  imagef([0 90],[360 -90],r);
  plotcont
  axis([280 350 50 88])
  caxis(gca,[0 prctile(r(:),95)]);
  colorbar
  longticks(gca,2)
  t=title('Spatial Variance of noise');
  
elseif strcmp(Clmlmp,'demo4')
  % Lets make a figure of the covariance of one point 
  % centered over Greenland
  L=60;
  [Clmlmp]=plmresid2cov([],L);
  [Covrr,Covrr_scaled,Covrr_shaped,rnot]=Clmlmp2Crrp(Clmlmp,[315 18],1);
  figure
  imagef([0 90],[360 -90],Covrr_shaped{1});
  plotcont
  % Put a white dot over out location of choice
  hold on
  h=plot(rnot(1,1),90-rnot(1,2),'wo');
  hold off
  set(h(1),'MarkerFaceColor','w')
  axis([0 360 -90 90])
  %caxis(gca,[0 1.5e-20]);
  colorbar
  longticks(gca,2)
  t=title('Spatial Covariance of noise');
elseif strcmp(Clmlmp,'demo5')
  % Same as demo4 but only spectral diagonals 
  L=60;
  [Clmlmp]=plmresid2cov([],L);
  [Covrr,Covrr_scaled,Covrr_shaped,rnot]=Clmlmp2Crrp(diag(diag(Clmlmp)),[315 18],1);
  figure
  imagef([0 90],[360 -90],Covrr_shaped{1});
  plotcont
  % Put a white dot over out location of choice
  hold on
  h=plot(rnot(1,1),90-rnot(1,2),'wo');
  hold off
  set(h(1),'MarkerFaceColor','w')
  axis([0 360 -90 90])
  %caxis(gca,[0 1.5e-20]);
  colorbar
  longticks(gca,2)
  t=title('Spatial Covariance of noise');

end % end if isstr

    
