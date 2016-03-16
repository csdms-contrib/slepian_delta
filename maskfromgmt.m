function varargout=maskfromgmt(XY,grdres,grdarea)
% [mymask]=MASKFROMGMT(XY,grdres,grdarea)
%
% This program takes a grid and a set of polygon regions, and uses the 
% program GRDMASK from the Generic Mapping Tools suite to generate a mask
% of your polygon areas.  This functions similar to Matlab's inpolygon
% routine but is much faster for very large datasets.
%
%
% INPUT:
% 
% XY        The polygons that you want to make a mask for.  These can be
%           multiple polygons separated by NaNs
% grdres    The resolution of the grid that you want to make a mask on. This 
%           should be string of the degree resolution [default: '1']
% grdarea   If you know the area where your polygons occur, you can save
%           yourself a lot of time by limiting your calculation to just
%           that area. This should be a string to be used in the GMT -R
%           parameter.
%
% OUTPUT:
%
% mymask        The mask from GMT as a N by 3 matrix, where the columns are
%               longitude, lattitude, and value. 1 is in the polygon, 0 is
%               outside.
% reshapedmask  Same mask as above, but reshaped into a grid for the
%               appropriate area
%     
% NOTE: This program expects version 5.X of GMT to be installed.
%
% SEE ALSO: GEOBOXCAP
%
% EXAMPLE:
% maskfromgmt('demo')
%
% Last modified by charig-at-princeton.edu, 03/16/2016

if ~isstr(XY)
    
    % Defaults
    defval('varlist',[]);
    defval('allsame',1);
    %defval('grdarea','260/300/74/84');
    defval('grdarea','0/360/-90/90');
    defval('grdres','1');

    %if ~ischar(grdres); grdres = num2str(grdres); end
    
    % Check that our polygons are the standard 0 to 360
    XY(:,1)=XY(:,1)-360*[XY(:,1)>360];

    % First we export our polygons to a text file
    fnpl1='mytemppolygons.dat';
    fp1 = fopen(fnpl1,'wt');
    fprintf(fp1,'%.4e %.4e \n',XY');
    fclose(fp1);

    % Replace the NaNs in that file with the GMT character for multiple polygons
    % For some reason we cannot do this in the same file
    system('perl -pe ''s/NaN NaN/>/g;'' mytemppolygons.dat > mytemppolygons2.dat');
    system('mv -f mytemppolygons2.dat mytemppolygons.dat');

    % Now we call gmt
    [a,b] = system(['gmt grdmask mytemppolygons.dat -Gmynewmask.grd -I' grdres ' -R' grdarea ' -V']);

    % Convert the grdfile back to a text file
    [a,b] = system(['grd2xyz mynewmask.grd > newmasktext.dat']);

    % Read this text file back into Matlab
    mymask = load('newmasktext.dat');

    % Clean up the files
    system('rm -f mynewmask.grd newmasktext.dat mytemppolygons.dat');

    % Reshape it
    [a,b,c,d] = strread(grdarea,'%n%n%n%n','delimiter','/');
    reshapedmask = reshape(mymask(:,3),length([a:str2num(grdres):b]),length([c:str2num(grdres):d]));
    
    % Collect output
    varns={mymask,reshapedmask};
    % Provide output where requested
    varargout=varns(1:nargout);

elseif strcmp(XY,'demo')
    XY = greenland(1);
    [mymask] = maskfromgmt(XY,'0.5','280/355/55/85');
    r = reshape(mymask(:,3),length([280:0.5:355]),length([55:0.5:85]));
    clf
    imagesc(r');
    
end

    
   
    
    
