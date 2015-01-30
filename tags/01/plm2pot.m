function lmcosi=plm2pot(lmcosi,r,GM,a,wat,wit)
% lmcosip=PLM2POT(lmcosi,r,GM,a,wat,wit)
%
% Converts the real spherical harmonic coefficients expressing a geoPOTENTIAL
% field in a convention that has GM/r pulled out of the sum to the
% corresponding coefficients for a particular gravitational field such as the
% geopotential, the free-air gravity, gravity disturbance, or the geoid
% anomaly.  Currently these are referenced to the WGS84 ellipsoid.
%
% If the (0,0) term is missing, the function adds in the reference
% constant for that as well, and also supplies the zeros for (1,0) and
% (1,1) to make the coordinate system explicitly the center-of-mass.
%
% INPUT:
%
% lmcosi     [l m Ccos Csin] degrees, order, coefficients, as in PLM2XYZ
%            Units are properly those of potential, J/kg or m^2/s^2
% r          Requested radius, in m [default: a]
% GM         Gravitational coefficient X mass reference [default: GM_EGM96]
% a          Reference radius [default: a_EGM96]
% wat        0 geoid anomaly with respect to 'wit' [m]
%            1 gravitational potential anomaly wrt 'wit' [J/kg]
%            2 free-air gravity anomaly wrt 'wit' [m/s^2]
%            3 gravity disturbance wrt 'wit' [m/s^2]
%            4 surface mass density [kg/m^2] from Wahr et al. (1998)
%              this assumes you have geoid coefficients
%            5 geoid [m] from surface mass density [kg/m^2]
%              (reverse of #4)
% wit        'nothing' applies to cases where no reference is subtracted
%            'WGS84' if your reference is WGS84 (even-degree zonals) [default]
%            'spherical' only ever takes out the (0,0) coefficients
%
% OUPUT:
%
% lmcosip    [l m Ccos Csin] degrees, order, coefficients of output
%            in the standard units of what's requested [J/kg, m/s^2, m]
%
% EXAMPLE:
%
% plm2pot('demo1') % Make field of some surface mass, convert to
% geoid, then convert back, making a plot of all 3.  Panel 1 and 3 should
% be the same.  Panel 2 should be more smoothed.
%
% Last modified by charig-at-princeton.edu, 02/23/2012
% Last modified by kwlewis-at-princeton.edu, 02/23/2012
% Last modified by fjsimons-at-alum.mit.edu, 05/24/2013

% Should write something that involves the geoid kernels and compares the
% dynamic ones with the load-Love ones. 

% See Research Book V page 84, Essentials 2001 and notes with egm96_?
% See Wahr et al. 1998, for surface density

if ~isstr(lmcosi)
    
  defval('GM',fralmanac('GM_EGM96','Earth'))
  defval('a',fralmanac('a_EGM96','Earth'))
  defval('r',a)
  defval('wat',1)
  defval('wit','WGS84')
  
  % Note this is always in here, don't confuse when using PLM2XYZ
  GMr=GM/r;
  
  % Blakely Eq (7.2)
  el=lmcosi(:,1);
  arl=(a/r).^el;
  
  % Published elastic Love numbers from Wahr et al. (1998).  
  % Perhaps worth updating these for a more realistic, but
  % still radial, mantle viscosity?
  % Interpolate the elastic Love numbers to the bandwidth you need
  [lkli,lkl]=lovenums('Wahr',el);
    
  switch wat
   case 0 % Geoid (anomaly) (undulation)
    fact=a;   
   case 1 % Potential (anomaly) - Blakely (7.2)
    fact=GMr*[arl arl];
    disp('Calculating gravitational potential')
   case 2 % Free air gravity anomaly, NOT the gravity disturbance!!!
    fact=GMr*[arl arl].*[el-1 el-1]/r;
    disp('Calculating free-air gravity anomalies')
   case 3 % Gravity disturbance, NOT the free air anomaly!!!
    fact=GMr*[arl arl].*[el+1 el+1]/r;
    disp('Calculating gravity disturbance')
   case 4
    % Surface mass density from geoid anomaly coefficients
    rhoave=5517; 
    disp(sprintf('Average Earth density used is %i kg/m^3',rhoave))
    fact=repmat(rhoave/3*(2*el(:)+1)./(1+lkli(:,2)),1,2);
    wit='nothing';
   case 5  
    % Reverse of 4, surface mass density to geoid anomaly
    rhoave=5517; 
    disp(sprintf('Average Earth density used is %i kg/m^3',rhoave))
    fact=repmat(3/rhoave*((1+lkli(:,2))./(2*el(:)+1)),1,2);
    wit='nothing';
  end
  
  switch wit
    case 'WGS84'
    % Make relative to the WGS 84 ellipsoid so calculate it geopotential first
    % Note that these are just [co l] instead of [l m co si]
    barC2n=grs([],[],[],[],5);
    % Supply the correct even zonal harmonic from the GRS earth model
    for l=[barC2n(:,2)]'
      were=addmup(l-1)-addmup(lmcosi(1,1)-1)+1;
      lmcosi(were,3)=lmcosi(were,3)-barC2n(l/2,1);
      disp(sprintf('Adjusting degree %2i with degree %2i',lmcosi(were,1),barC2n(l/2,2)))
    end
    disp('Referenced to the WGS-84 ellipsoid')
  end
  
  % Complete the potential coefficients
  if lmcosi(1)==2
    % Mean and center of mass were zero - supply them pro forma
    lmcosi=[0 0 GMr 0 repmat(0,1,size(lmcosi,2)-4); % Degree 0
	    1 0 0 0 repmat(0,1,size(lmcosi,2)-4); % Degree 1 order 0
	    1 1 0 0 repmat(0,1,size(lmcosi,2)-4); % Degree 1 order 1
	    lmcosi]; % Then the rest
  elseif lmcosi(1)==1 && lmcosi(1,2)==0
    lmcosi=[0 0 GMr 0 repmat(0,1,size(lmcosi,2)-4); % Degree 0
	    lmcosi]; % Then the rest
  end

  switch wit
   case {'WGS84','spherical'}
    % Take out the (0,0) reference again, after all
    lmcosi(1,3)=0;
  end
  
  % Factor in the factor for everything else after this
  lmcosi(:,3:4)=fact.*lmcosi(:,3:4);
  
elseif strcmp(lmcosi,'demo1')
  % Make some surface mass in a region easy to plot
  [V,C,dels,dems,XY]=localization(18,'samerica');
  C{1}=C{1}*2.94; % Scale by 2.94 so the peak is roughly 50 kg/m^2 (like avg Greenland Melt)
  % Say the best function is our surface mass and convert it to geoid.
  [lmcosiG]=plm2pot([dels dems C{1}],[],[],[],5);
  % Now put it back in surface density
  [lmcosiSD]=plm2pot(lmcosiG,[],[],[],4);
  % Now make a figure of these 3 things
  ah=krijetem(subnum(3,1));
  fig2print(gcf,'tall')
  % Plot the spatial data in Mollweide projection
  axes(ah(1))
  [data,ch,ph]=plotplm([dels dems C{1}],[],[],4,1);
  kelicol
  t=title('Slepian fn from LOCALIZATION representing our surface density [kg/m^2]');
  %caxis([-13 1])
  cb=colorbar('hor');
  shrink(cb,1.5,2)
  movev(cb,-0.05)
  axes(cb)
  set(get(cb,'xlabel'),'string','Surface density [kg/m^2]');
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  axes(ah(2))
  % Plot in mm
  lmcosiG(:,3:4)=lmcosiG(:,3:4)*1000;
  [data,ch,ph]=plotplm(lmcosiG,[],[],4,1);
  kelicol
  t=title('Panel 1 represented as Geoid [milimeters]');
  %caxis([-0.5 0.05])
  cb=colorbar('hor');
  shrink(cb,1.5,2)
  movev(cb,-0.05)
  axes(cb)
  set(get(cb,'xlabel'),'string','Geoid height variation [milimeters]');
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  axes(ah(3))
  [data,ch,ph]=plotplm(lmcosiSD,[],[],4,1);
  kelicol
  t=title('Panel 2 converted back to surface density [kg/m^2]');
  %caxis([-13 1])
  cb=colorbar('hor');
  shrink(cb,1.5,2)
  movev(cb,-0.05)
  axes(cb)
  set(get(cb,'xlabel'),'string','Surface density [kg/m^2]');
end

