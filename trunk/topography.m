function topography
% TOPOGRAPHY
%
% Illustrates Earth, Moon, Mars topography by plotting it in space as
% well as by calculating its spherical harmonic spectrum.
%
% SEE ALSO: FRS149_1, FRS149_3
%
% Last modified by fjsimons-at-alum.mit.edu, Feb 12th, 2004

defval('wat','XYZ')

P{1}=fralmanac('GTM3AR',wat);     % Earth SHM or XYZ
U{1}=fralmanac('GTM3AR','SHM');   % Earth SHM or XYZ

%P{2}=fralmanac('GLTM2B','SHM');   % Moon
%U{2}=P{2};
%P{3}=fralmanac('GTM090','SHM');   % Mars
%U{3}=P{3};

if length(P)>=3
  % For Mars, take out degree two to see more
  P{3}(4,3)=0;
end
N={'Earth' 'Moon' 'Mars'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clf
[ah,ha]=krijetem(subnum(length(P),2));

for index=1:length(P)
  % Plot field
  axes(ha(index))
  [data,ch{index},ph{index}]=plotplm(P{index},[],[],1); 
  tl(index)=title(N{index});
  cb(index)=colorbar('hor');
 
  axes(ha(index+length(P)))
  % Plot spectrum
  splo{index}=plotplm(U{index},[],[],3);
  set(splo{1}(2),'LineW',2)
  yl(index)=ylabel('Spectral Density');
  xl(index)=xlabel('Degree');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Delete continents and plates from Mars and Moon
if length(P)>1
  delete([ch{2}{1} ch{3}{1} ph{2} ph{3}])
  shrink(cb,1.5,1)
  shrink(ha(length(P)+1:end),1,1.3)
end

fig2print(gcf,'landscape')
id
figdisp
disp('Use option -painters')


 
