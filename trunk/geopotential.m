function varargout=geopotential(wat)
% data=GEOPOTENTIAL(wat)
%
% Illustrates Earth, Moon, Mars gravity
%
% INPUT:
%
% wat    A vector with, for Earth, Moon and Mars:
%            1 gravitational potential
%            2 free-air gravity anomaly (default)
%            3 geoid anomaly
%
% OUTPUT:
%
% data   The data being plotted
%
% SEE ALSO:
%
% FRS149_1, FRS149_3, PLM2POT, EGM96_X
%
% Last modified by fjsimons-at-alum.mit.edu, 06/19/2008

defval('wat',2)

P{1}=plm2pot(fralmanac('EGM96','SHM'),...
	     fralmanac('a_EGM96','Earth'),...
	     fralmanac('GM_EGM96','Earth'),...
	     fralmanac('a_EGM96','Earth'),...
	     wat(1));
if length(wat)>1
  P{2}=plm2pot(fralmanac('GLGM2','SHM'),...
	       fralmanac('a_GLGM2','Moon'),...
	       fralmanac('GM_GLGM2','Moon'),...
	       fralmanac('a_GLGM2','Moon'),...
	       wat(2));
  if length(wat)>2
    P{3}=plm2pot(fralmanac('JGM85H02','SHM'),...
		 fralmanac('a_JGM85H02','Mars'),...
		 fralmanac('GM_JGM85H02','Mars'),...
		 fralmanac('a_JGM85H02','Mars'),...
		 wat(3));
  end
end

N={'Earth','Moon','Mars'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clf
[ah,ha]=krijetem(subnum(length(P),2));

for index=1:length(P)
  axes(ha(index))
  L=addmup(length(P{index}),'r');
  % This is the Nyquist pixel size warranted by the maximum degree
  degN=180/sqrt(L*(L+1));
  % Makes a plot of the field generated from the coefficients
  if wat(index)==2
    % Convert to mgal
    P{index}(:,3:4)=P{index}(:,3:4)*1e5;
  end
  [data{index},ch{index},ph{index}]=plotplm(P{index},[],[],1,degN); 
  if wat(index)==2
    caxis([-100 100])
  end
  tl(index)=title(N{index});
  cb(index)=colorbar('SouthOutside');
  if getpos(cb(index),4)<0
    set(cb(index),'position',[getpos(cb(index),1:3) 0.03])
  end
  axes(ha(index+length(P)))
  % Makes a plot of the spectrum of this field: total energy per degree
  splo{index}=plotplm(P{index},[],[],3);
  yl(index)=ylabel('Spectral Density');
  pxl(index)=xlabel('Degree');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Delete continents and plates from Mars and Moon
if length(P)>1
  delete([ch{2}{1} ph{2}])
  if length(P)>2
    delete([ch{3}{1} ph{3}])
    shrink(ha(length(P)+1:end),1,1.3)
    shrink(cb,1.5,1)
  end
  fig2print(gcf,'tall')
else
  shrink(ah(2),1,2.5)
  % Rather put the bottom of the colorbar where the bottom of the rhs is
  % movev(cb,.2)
  hh=getpos(cb);
  set(cb,'Position',[hh(1) getpos(ah(2),2)-0.025 hh(3:4)])
  fig2print(gcf,'portrait')
end

figdisp

varns={data};
varargout=varns(1:nargout);

