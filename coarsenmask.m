function varargout=coarsenmask(degres,c11cmn,fine,method,option1,option2)
% [mymask]=COARSENMASK(degres,c11cmn,fine,method,option1,option2)
% [mymask]=COARSENMASK(coarselon,coarselat,fine,method,option1,option2)
%
% This program takes a mask, such as a land mask, and coarsens the
% resolution by a few methods.
%
%
% INPUT:
% 
% degres     Longitude/ latitude spacing, in degrees [default: Nyquist] OR
%             "lat": a column vector with latitudes [degrees]
% c11cmn     Corner nodes of lon/lat grid [default: 0 90 360 -90] OR
%             "lon": a column vector with longitudes [degrees]
% fine       The finer resolution mask you want to coarsen.  This can be
%            either a full matrix, or perhaps some vectors of coordinates.
%            For example you can send it just the parts that have mask=1 to
%            save time.  In this case you also need to give a dlon depending
%            on what option you choose.
% method     [1] Any part of the old mask within the box of the new mask 
%                is transferred [default]
%            [2] Parts of the old mask are area weighted in the new mask
%            [3] A certain % (by area) of the old mask has to be 1 for the
%                new mask to be 1
% option1    The resolution of the fine mask that you are giving the 
%            program.  Lon and Lat resolution is assumed to be the same.
%            Used in method 2 and 3.
% option2    The threshold for area weighting.  0.5 [default] means that
%            50% of the coarse box by area must contain "1" in the fine 
%            mask, otherwise it is set to zero. Used in method 3.  
%
% OUTPUT:
%
% mymask    The new coarsemask as a m-by-n matrix, similar to the output of
%           plm2xyz.
% lonlon    The Cartesion x-coordinates of the geographical [lon lat] grid
% latlat    The Cartesion y-coordinates of the geographical [lon lat] grid
%     
%
% EXAMPLE:
% coarsenmask('demo')
%
% Last modified by charig-at-princeton.edu, 03/16/2016

defval('method',1)
defval('degres',2);
defval('c11cmn',[0 90 360 -90]);

if ~isstr(degres) % Not a demo
    
    % Do we have a square and resolution or vectors for the lon/lat?
    if length(degres)==1 && length(c11cmn)==4
        % It's a grid
        c11 = c11cmn(1:2);
        cmn = c11cmn(3:4);
        coarselon = [c11cmn(1):degres:c11cmn(3)];
        coarselat = [c11cmn(4):degres:c11cmn(2)];
    elseif length(degres)>1 && length(c11cmn)~=4
        % It must be vectors for Lon and Lat, such as from geoboxcap
        % Check that coarselat is descending
        coarselon = degres;
        coarselat = c11cmn;
        coarselat = sort(coarselat,'descend');
        c11 = [min(coarselon) max(coarselat)];
        cmn = [max(coarselon) min(coarselat)];
    end

    
    % What is fine?
    if iscell(fine)
        % You are passing some output from GEOBOXCAP maybe
        % in which case the lon and lat are row vectors
        finemask = fine{1};
        finemasklon = gamini(fine{2},length(fine{3}))'; % one for each lon, then flip
        finemasklat = repmat(fine{3}',length(fine{2}),1);
    elseif size(fine,2)==3
        % You passed it coordinates, such as from MASKFROMGMT
        finemasklon = fine(:,1);
        finemasklat = fine(:,2);
        finemask = fine(:,3);
    else
        error('Something wrong with FINE?')
    end
    

    m = length(coarselat);
    n = length(coarselon);
    % The coarsegrid resolution
    lont=indeks(diff(coarselon),1); latt=abs(indeks(diff(coarselat),1));
    % Allocate the new coarsemask
    [lonlon,latlat]=meshgrid(coarselon,coarselat);
    latlat = flipud(latlat);
    coarsemask = zeros(size(lonlon));

    % Only grab the parts that are 1 in the old mask.
    finemasklon = finemasklon(finemask==1);
    finemasklat = finemasklat(finemask==1);
    
    % Now figure out where the fine mask is in relation to the coarse mask.
    % We can use COR2IND, however we need to account for the fact that this
    % function is meant to work with grids which have smaller sizes than
    % their underlying matrices
    %c11 = [min(coarselon)-0.5 max(coarselat)+0.5];
    %cmn = [max(coarselon)+0.5 min(coarselat)-0.5];
    [ind,colnr,rownr] = cor2ind(finemasklon,finemasklat,...
        [c11(1)-lont/2 c11(2)+latt/2],[cmn(1)+lont/2 cmn(2)-latt/2],m,n);
    
    % Do the method we want
    if method==1
        % Any part of the old mask within the box of the new mask 
        % is transferred [default]
        coarsemask(ind) = 1;
        
    elseif method==2
        % Parts of the old mask are area weighted in the new mask
        defval('option1','indeks(diff(finemasklon),1)') % Try this
        % Get the areas represented by the coarse grid points
        c11_c = [lonlon(:)-lont/2 latlat(:)+latt/2];
        cmn_c = [lonlon(:)+lont/2 latlat(:)-latt/2];
        coarseareas = spharea(c11_c,cmn_c);
        
        % Get the areas represented by the fine grid points
        dlon=option1; dlat=dlon;
        c11_f = [finemasklon-dlon/2 finemasklat+dlat/2];
        cmn_f = [finemasklon+dlon/2 finemasklat-dlat/2];
        fineareas = spharea(c11_f,cmn_f);
        
        % Calculate how much area the fine mask takes up in each coarse box
        fineareasums = accumarray(ind,fineareas);
        
        % Divide this by the total coarse area to get the area weight
        coarsemask(ind) = fineareasums(ind)./coarseareas(ind);
        
    elseif method==3
        % A certain % (by area) of the old mask has to be 1 for the new mask to be 1
        defval('option1','indeks(diff(finemasklon),1)') % Try this
        % Do the same thing as method 2, but just round at the end
        
        % Get the areas represented by the coarse grid points
        c11_c = [lonlon(:)-lont/2 latlat(:)+latt/2];
        cmn_c = [lonlon(:)+lont/2 latlat(:)-latt/2];
        coarseareas = spharea(c11_c,cmn_c);
        % Get the areas represented by the fine grid points
        dlon=option1; dlat=dlon;
        c11_f = [finemasklon-dlon/2 finemasklat+dlat/2];
        cmn_f = [finemasklon+dlon/2 finemasklat-dlat/2];
        fineareas = spharea(c11_f,cmn_f);
        % Calculate how much area the fine mask takes up in each coarse box
        fineareasums = accumarray(ind,fineareas);
        % Divide this by the total coarse area to get the area weight
        coarsemask(ind) = fineareasums(ind)./coarseareas(ind);
    
        % Threshold!
        defval('option2','0.5')
        coarsemask(coarsemask<option2) = 0;
    
    end
    
    % Collect output
    varns={coarsemask,lonlon,latlat};
    % Provide output where requested
    varargout=varns(1:nargout);

elseif strcmp(degres,'demo')

    XY = namerica();
    fineres = 0.5;
    [mymask] = maskfromgmt(XY,'.5');
    coarselon = [0:1:360];
    coarselat = [-90:1:90];
    
    figure
    clf
    fig2print(gcf,'portrait')
    ah1=krijetem(subnum(3,1));
    
    axes(ah1(1));
    temp1 = coarsenmask(coarselon,coarselat,mymask,1);
    imagesc(temp1)
    colorbar
    axis equal
    axis([190 310 10 85])
    title('Method 1: Any part goes to new mask');
%     XY = namerica();
%     hold on
%     plot(XY(:,1),90-XY(:,2),'w')
    
    axes(ah1(2));
    temp1 = coarsenmask(coarselon,coarselat,mymask,2,fineres);
    imagesc(temp1)
    colorbar
    axis equal
    axis([190 310 10 85])
    title('Method 2: Area weighting');
    
    
    axes(ah1(3));
    temp1 = coarsenmask(coarselon,coarselat,mymask,3,fineres,0.7);
    imagesc(temp1)
    colorbar
    axis equal
    axis([190 310 10 85])
    title('Method 3: 70% threshold area weight');
    
    
end

    
   
    
    
