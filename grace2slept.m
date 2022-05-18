function varargout=grace2slept(Dataproduct,TH,XY_buffer,Lwindow,phi,theta,omega,J,units,forcenew)
% [slepcoffs,calerrors,thedates,TH,G,CC,V]
%                 =GRACE2SLEPT(Dataproduct,TH,XY_buffer,Lwindow,phi,theta,omega,J,units,forcenew)
%
% This program takes GRACE/GRACE-FO gravimetry data created by GRACE2PLM
% and projects this data into the requested Slepian basis.
%
%
% INPUT:
% 
% Dataproduct   This is a cell array with three parts: a) the data center,
%                 b) the release level, and c) the dataproduct bandwidth. 
%                 [default: {'CSR', 'RL06', 60} which is Release level 06 
%                 data from the data center]
%                 at the Center for Space Research
%                 Another Example: {'GFZ', 'RL05', 96} which is Release 
%                 level 05 data from the data center at the 
%                 GeoForschungsZentrum Potsdam
% TH         Radius of the concentration region (degrees) OR
%              'england', 'eurasia',  'namerica', 'australia', 'greenland'
%              'africa', 'samerica', 'amazon', 'orinoco', in which case 
%              you must have phi,theta,omega all equal to zero OR
%              [lon lat] an ordered list defining a closed curve [degrees]
% XY_buffer  Distance in degrees that the region outline will be enlarged
%              by BUFFERM [default: 0]
% Lwindow    Bandwidth of the window [default: bandwidth of the data], 
%              or bandpass (two degrees)
% phi        Longitude of the center of the tapers (degrees)
% theta      Colatitude of the center of the tapers (degrees)
% omega      Anticlockwise azimuthal rotation of the tapers (degrees)
% J          Number of largest eigenfunctions in which to expand 
%             [default: all of them]  Can give the string 'N' if you want
%             the rounded Shannon number.
% units      'POT' or 'SD' for whether you want geopotential or surface
%            mass density
% forcenew    Wether or not you want to force new generation of a save file
%              (1) or just use the one we already have (0) [default].
%
% OUTPUT:
% 
% Returns these variables and saves the first three in a .mat 
% file (kernels are not saved):
%  slepcoffs    The expansion coefficients of the geopotential (or surface 
%                   density) into the Slepian basispotential Slepian coefficients 
%                   [nmonths x addmoff(Ldata)]
%  slepcalerrors  The expansion coefficients of the calibrated errors 
%                   into the Slepian basis calibrated errors 
%                   [nmonths x addmoff(Ldata)]
%  thedates     Time stamps in Matlab time
%  TH           The region back to you.  If there was buffering, this will
%                 be a XY array of coordinates, which you can use with 
%                 SPHAREA to get the Shannon number.
%  G            The unitary matrix of localization coefficients
%  CC           A cell array with cosine/sine coefficients eigenfunctions
%  V            The eigenvalues in this ordering
%  N            The Shannon number
%
% SEE ALSO: PLM2SLEP
%
%
% Last modified by charig-at-princeton.edu, 05/18/2022
% Last modified by fjsimons-at-alum.mit.edu, 06/26/2012

% Determine parameters and set defaults
defval('Dataproduct',{'CSR', 'RL06', 60})
defval('TH','greenland')
defval('Lwindow',18)
defval('phi',0)
defval('theta',0)
defval('omega',0)
defval('forcenew',0)
defval('XY_buffer',0)
defval('pars',10);
defval('units','POT')
defval('inout','out')
Pcenter = Dataproduct{1};
Rlevel = Dataproduct{2};
Ldata = Dataproduct{3};
dataproductstring = [Pcenter Rlevel num2str(Ldata)];

% Figure out if it's lowpass or bandpass
lp=length(Lwindow)==1;
bp=length(Lwindow)==2;
maxL=max(Lwindow);
% The spherical harmonic dimension
ldim=(Lwindow(2-lp)+1)^2-bp*Lwindow(1)^2;
defval('J',ldim)

% Where the original data files are kept
defval('ddir1',fullfile(getenv('IFILES'),'GRACE','Originals',Rlevel,Pcenter));
% Where you would like to save the new .mat file
defval('ddir2',fullfile(getenv('IFILES'),'GRACE','SlepianExpansions'));

% Get the remaining file names
if ~isstr(TH) && length(TH)==1 % POLAR CAPS
    % Check if you want the Shannon number of eigenfunctions
    if strcmp(J,'N')
       J = round((Lwindow+1)^2*spharea(TH));
    end
    if lp
        fnpl=sprintf('%s/grace2slept-%s-CAP-%i-%i-%i-%s.mat',...
            ddir2,dataproductstring,TH,Lwindow,J,units);
    elseif bp
        fnpl=sprintf('%s/grace2sleptbl-%s-CAP-%i-%i-%i-%i-%s.mat',...
            ddir2,dataproductstring,TH,Lwindow(1),Lwindow(2),J,units);
    else
        error('The degree range is either one or two numbers')       
    end 
        
else % GEOGRAPHICAL REGIONS and XY REGIONS
    % Check if you want the Shannon number of eigenfunctions
    if strcmp(J,'N')
       J = round((Lwindow+1)^2*spharea(TH));
    end
    if isstr(TH) % Geographic
        % Here, TH gets passed to glmalpha, and glmalpha will interpret
        % either the cell of the region
        if XY_buffer ~= 0
            TH = {TH XY_buffer};
            h = [TH{1} num2str(XY_buffer)];
        else
            h=TH;
        end
    else % Closed coordinates (make a hash)
        h=hash(TH,'sha1');
    end
    % The name of the save file
    if lp
        fnpl=sprintf('%s/grace2slept-%s-%s-%i-%i-%s.mat',...
            ddir2,dataproductstring,h,Lwindow,J,units);
    elseif bp
        fnpl=sprintf('%s/grace2sleptbl-%s-%s-%i-%i-%i-%s.mat',...
            ddir2,dataproductstring,h,Lwindow(1),Lwindow(2),J,units);
    else
        error('The degree range is either one or two numbers')
    end

end

% GET THE SLEPIAN BASIS WE WANT
[~,~,~,lmcosiW,~,~,~,~,~,ronmW]=addmon(maxL);
% NOTE: a) Could have used PLM2SLEP but since we need the same operation for many
%         months, slightly better to load the Slepian basis once, and repeatedly
%         multiply by G.  Mostly copied from PLM2SLEP.
%       b) The kernel and eigenfunctions are large to save for large Lwindow.  
%          So instead just load them from either GLMALPHA or GLMALPHAPTO.
%       c) GLMALPHA now handles bandpass

% If it is the standard North-Polar cap or a geographic region, it's easy
if phi==0 && theta==0 && omega==0    
   % Get the Slepian basis; definitely not block-sorted as for the rotated
   % versions this will make no sense at all anymore
   [G,V,EL,EM,N,GM2AL,MTAP,IMTAP]=glmalpha(TH,Lwindow,[],0,[],[],J);
else
   % Need to get a complete GLMALPHA but for the rotated basis
   % Definitely, "single-order" has lost its meaning here, but the MTAP
   % will still identify what the order of the unrotated original was
   [G,V,EL,EM,N,GM2AL,MTAP,IMTAP]=glmalphapto(TH,Lwindow,phi,theta,omega);
   % Since GLMALPHAPTO currently has no option to limit a basis to J, do it here
   G = G(:,1:J);
end
% Sort by decreasing eigenvalue
[V,vi]=sort(V,'descend');
G=G(:,vi); if ~isnan(MTAP); MTAP=MTAP(vi); end
% If you don't do this, the eigenfunctions are ordered in the way
%   that they correspond to single-orders back when, unrotated, they
%   belonged to a polar cap, and the eigenvalues are sorted within
%   these blocks. This is useful for, e.g. SPIE2009_1 a la SDSNEEUW. 
% Collect the eigenvector output into a format that PLM2XYZ knows how to interpret
for j=1:size(G,2)
   % Create the blanks
   cosi=lmcosiW(:,3:4);
   % Stick in the coefficients of the 1st eigentaper
   cosi(ronmW)=G(:,j);
   % Construct the full matrix
   CC{j} = [lmcosiW(:,1:2) cosi]; 
end


% INITILIZATION COMPLETE

% If this expansion already exists, load it.  Otherwise, or if we force 
% it, make a new one (e.g. if you added extra months to the database).
if exist(fnpl,'file')==2 && forcenew==0
     load(fnpl)
     disp(sprintf('%s loaded by GRACE2SLEPT',fnpl))
else
    % Use GRACE2PLMT to get the GRACE data.  This way, if we have it saved,
    % there is no need to scan the month files again.  GRACE2PLMT takes 
    % care of the WGS84 adjustment, the degree 1 correction, and the C20 correction. 
    [potcoffs,cal_errors,thedates]=grace2plmt(Pcenter,Rlevel,Ldata,units,0);
    % *** Here I changed this. Run grace2plmt once to update your data, and
    % then when you call forcenew=1 from now on it will just update the
    % expansion
    
    % Initialize new coefficients
    nmonths=length(thedates);
    slepcoffs=nan(nmonths,J);
    slepcalerrors=nan(nmonths,J);

    % Limit everything to the window bandwidth
    potcoffsW = potcoffs(:,1:size(lmcosiW,1),1:4);
    %cal_errorsW = cal_errors(:,1:size(lmcosiW,1),1:4);

    % Loop over the months
    for index = 1:nmonths
        % Expand this month's POTENTIAL into the Slepian basis
        potcoffs_month=squeeze(potcoffsW(index,:,:));
        slepcoffs(index,:) = ...
            potcoffs_month(2*size(potcoffs_month,1)+ronmW(1:(maxL+1)^2))'*G;

        % Expand this month of CALIBRATED ERRORS into the Slepian basis
        %calerrors_month=squeeze(cal_errorsW(index,:,:));
        %slepcalerrors(index,:) = ...
        %    calerrors_month(2*size(calerrors_month,1)+ronmW(1:(maxL+1)^2))'*G;
    end
       

    % SAVE
    % Here we have "thedates" twice so that we don't break older code. But in
    % the future we will fix this so that we don't have cal errors
    % Don't save the kernel and eigenfunctions because we already have 
    % this info saved and can load from GLMALPHA
    % save(fnpl,'slepcoffs','calerrors','thedates','G','CC','V');
    save(fnpl,'slepcoffs','thedates','thedates');

end % End if we have a save file already

% Collect output
varns={slepcoffs,thedates,thedates,TH,G,CC,V,N};
varargout=varns(1:nargout);
