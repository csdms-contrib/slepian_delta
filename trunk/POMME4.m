function varargout=POMME4(L,actplot,degres)
% [d,lmcosip,degres]=POMME4(L,actplot,degres)
%
% Plots a lithospheric magnetic field model, POMME, which is complete
% from degree 1 and order 0 to degree and order 720. 
%
% INPUT:
%
% L          Truncation degrees: 
%               if NaN, just returns the map
%               if empty, default is [17 72]
% actplot    1 Actually plot this [default]
%            0 Just return the data
% degres     The degree resolution
%
% OUTPUT:
%
% d          The map being plotted
% lmcosip    The spherical harmonics being plotted
% degres     The degree resolution
%
% Makes a map of the magnetic field, and returns the data
% plotted if requested.
%
% Last modified by fjsimons-at-alum.mit.edu, 02/18/2011

defval('actplot',1)
defval('degres',1/4)
defval('L',[])

% Has it been expanded previously?
if isnan(L)
  filename=fullfile(getenv('IFILES'),...
		    'EARTHMODELS','POMME-4','POMME-4_BrnT.mat');
else
  defval('L',[17 72])
  filename=fullfile(getenv('IFILES'),...
		    'EARTHMODELS','POMME-4',...
		    sprintf('POMME-4_BrnT_%i_%i-%i.mat',L(1),L(2),degres));
end

if exist(filename,'file')~=2
  % Load the coefficients
  lmcosi=load(fullfile(getenv('IFILES'),...
		       'EARTHMODELS','POMME-4','pomme-4.2s-nosecular.cof'));
  if ~isempty(L)
    missl=addmup(lmcosi(1,1)-1);
    lmcosi=lmcosi(addmup(L(1)-1)-missl+1:addmup(L(2))-missl,:);
  end
  
  % Figure out if the dimensions are right
  lp=length(L)==1; bp=length(L)==2;
  ldim=(L(2-lp)+1)^2-bp*L(1)^2;
  difer(ldim-(2*length(lmcosi)-[L(2)-L(1)]-1))
  
  % Convert to radial-component magnetic field on the reference surface
  lmcosip=plm2mag(lmcosi);

  % Then expand (and plot)
  if actplot
    clf
    [d,ch,ph]=plotplm(lmcosip,[],[],4,degres);
  else
    d=plm2xyz(lmcosip,degres);
  end
  
  % Then save
  save(filename,'d','lmcosip')
else
  load(filename)
  disp(sprintf('Loading %s',filename))
  if actplot
    clf
    imagef(d); kelicol; plotcont; axis image; caxis(halverange(d))
    longticks(gca,2); hold off; cb=colorbar('hor'); 
    movev(cb,-0.1); longticks(cb,2); 
    axes(cb); hold off
    xlabel(['Lithospheric magnetic field (nT) for degrees' ...
	    ' l=16...720'])
    fig2print(gcf,'landscape')
  end
end

% Output if requested
vars={d,lmcosip,degres};
varargout=vars(1:nargout);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cols=halverange(data)
cols=[-max(abs(data(:))) max(abs(data(:)))]/2;
