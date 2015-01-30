function plmt2diff(month1,month2,Lf)
% PLMT2DIFF(month1,month2,Lf)
%
% Makes a map of the (filtered) free-air anomaly differences between two
% monthly indices in the vector plmt created by GRACE2PLMT. If the months
% are the same, just plot the single month already.
%
% Sumatra? 30,31
% Chile? 91,92
%
% Last modified by fjsimons-at-alum.mit.edu, 01/14/2013

defval('month1',82);
defval('month2',88);
% Filtering
defval('Lf',18);

clear plmt
defval('plmt',fullfile(getenv('IFILES'),'GRACE','CSR_alldata.mat'))
a=load(plmt);
% These are the potential coefficients
plmt=a.potcoffs(:,:,1:4);
% Now get the thedates for each of those months
thedates=a.thedates;
clear a

% Load the auxiliary parameters, really should do this from the header
a=fralmanac('a_EGM2008','Earth');
GM=fralmanac('GM_EGM2008','Earth');

% Pick two months, use DATEVEC, DATESTR, DATENUM to convert
thedates=thedates([month1 month2]);
plmt=plmt([month1 month2],:,:);

% Normalize to give the gravity at a as per readme.egm96 and Blakely
GMa=GM/a;

% The degrees
el=squeeze(plmt(1,:,1));

% This already wrt WGS84 (and the 20 was being measured by SLR)

% But now we take out the DC term
plmt(1,1,[3 4])=0;
plmt(2,1,[3 4])=0;

% The following factor gives us the anomalous "free-air"
fact=GMa.*[el-1]/a;

% First selected month, cosine coefficient
plmt(1,:,3)=plmt(1,:,3).*fact;
% First selected month, sine coefficient
plmt(1,:,4)=plmt(1,:,4).*fact;
% Second selected month, cosine coefficient
plmt(2,:,3)=plmt(2,:,3).*fact;
% Second selected month, sine coefficient
plmt(2,:,4)=plmt(2,:,4).*fact;

% Reshape to standard lmcosi format
plmt1=reshape(plmt(1,:,[1 2 3 4]),[],4);
plmt2=reshape(plmt(2,:,[1 2 3 4]),[],4);

% Full-resolution of EGM2008
L=1500;
degres=180/sqrt(L*(L+1));

% Expand to space after filtering
[r1,lon,lat]=plm2xyz(plmfilt(plmt1,Lf),degres);
[r2,lon,lat]=plm2xyz(plmfilt(plmt2,Lf),degres);

clf
ah=gca;

% Plot on the graph
if month1==month2
  r1=0;
  xs='free-air anomaly';
  modx=sprintf('in %s',datestr(datevec(thedates(2)),'mmmm yyyy'));
else
  xs='difference in free-air anomaly';
  modx=sprintf('between %s and %s',...
	       datestr(datevec(thedates(2)),'mmmm yyyy'),...
	       datestr(datevec(thedates(1)),'mmmm yyyy'));
end
[rdiff,c,ph]=plotplm(setnans(r2-r1));

% Reasonable color scale
% Times 5 for Sumatra
% Times 5000 for single month
caxis([-1 1]*1e2*1e-9*1.25)

kelicol
tt=kelicol;
tt(81,:)=[1 1 1];
tt(82,:)=[1 1 1];
colormap(tt)

cb=colorbar('hor');
axes(cb)
xlabel(sprintf(...
    '%s filtered to L = %i [m/s^2]\n %s %s',...
    xs,Lf,modx))
set(cb,'xaxisl','t')

shrink(cb,2,2)
axes(cb)
longticks(cb,2)
movev(cb,-.05*2)

set(ah,'camerav',6)
movev(cb,-.05*2)
movev([ah cb],.1)

moveh([ah cb],-.015)
%print -djpeg -r300 /home/fjsimons/EPS/alumnitea2_300
%print -djpeg -r600 /home/fjsimons/EPS/alumnitea2_600
%print -djpeg -r300 /home/fjsimons/EPS/climatechange
