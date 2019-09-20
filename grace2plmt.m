function varargout=grace2plmt(Pcenter,Rlevel,units,forcenew)
% [potcoffs,cal_errors,thedates]=GRACE2PLMT(Pcenter,Rlevel,units,forcenew)
%
% This program reads in the Level-2 GRACE geoid products from either the CSR or
% GFZ data centers, does some processing, and saves them as a plmt matrix
% in a .mat file.  In particular, the coefficients are reordered to our
% prefered lmcosi format, they are referenced to the WGS84 ellipsoid, 
% the C2,0 coefficients are replaced with more accurate measurements from
% satellite laser ranging, and the degree one coefficients are 
% substituted with those from Swenson et al. (2008).  You have the option 
% of leaving them as geopotential 
% or converting them to surface mass density using the method of 
% Wahr et al. 1998, based on Love numbers (see PLM2POT).
%
% INPUT:
% 
% Pcenter     'CSR' data center at the Center for Space Research
%             'GFZ' data center at the GeoForschungsZentrum Potsdam
% Rlevel      The release level of the solution you want.  
%              Either 'RL04','RL05', or 'RL06'
% units       'POT' or 'SD' for whether you want geopotential or surface
%               mass density
% forcenew    Whether or not you want to force new generation of a save file
%              (1) or just use the one we already have (0) [default].
%
% OUTPUT:
% 
% Returns these variables and saves them in a .mat file:
%    potcoffs       potential coefficients [nmonths x addmup(Ldata) x 6]
%                    these could also be in surface mass density
%    cal_errors     calibrated errors [nmonths x addmup(Ldata) x 4]
%    thedates       time stamps in Matlab time
%
% NOTE:

%	SLR data available from the GRACE Tellus website for RL06:
%	ftp://podaac-ftp.jpl.nasa.gov/allData/grace/L2/CSR/RL06
%	 
%	The RL05 solutions from CSR do not have any standard deviations 
%    given, formal or calibrated.  These values will be reported as 0.
%
%   SLR data available from the GRACE Tellus website:
%   http://grace.jpl.nasa.gov/data/J2/ notably 
%   ftp://ftp.csr.utexas.edu/pub/slr/degree_2/C20_RL04_2010_12.txt
%  The header was removed and the file renamed for easy use.
%  Updated files keep getting posted in the same location.  
%
%  SLR data available from the GRACE Tellus website:
%   http://grace.jpl.nasa.gov/data/degree1/ notably 
%   ftp://podaac.jpl.nasa.gov/allData/tellus/L2/degree_1/
%  The header was removed and the file renamed for easy use.
%  Updated files keep getting posted in the same location.  
%
% EXAMPLE: to make a new save file when you have added more months
% [potcoffs,cal_errors,thedates]=grace2plmt('CSR','RL05','SD',1);
%
%
% Last modified by mlubeck-at-email.arizona.edu, 09/12/2019
% Last modified by charig-at-email.arizona.edu, 03/16/2016
% Last modified by fjsimons-at-alum.mit.edu, 05/17/2011

% Determine parameters and set defaults
defval('Pcenter','CSR')
defval('Rlevel','RL06')
defval('units','SD')
defval('forcenew',1)

% Where the original data files are kept
defval('ddir1',fullfile(getenv('IFILES'),'GRACE','Originals',Rlevel,Pcenter));
 

% Where you would like to save the new .mat file
defval('ddir2',fullfile(getenv('IFILES'),'GRACE'));
% And the name of that save file
if strcmp(units,'SD')
    fnpl=sprintf('%s/%s_%s_alldata_%s.mat',ddir2,Pcenter,Rlevel,units);
elseif strcmp(units,'POT')
    fnpl=sprintf('%s/%s_%s_alldata.mat',ddir2,Pcenter,Rlevel);
else
    error ('The unit input is not valid.')
end

% If this file already exists, load it.  Otherwise, or if we force it, make
% a new one (e.g. you added extra months to the database).
if exist(fnpl,'file')==2 && forcenew==0
     load(fnpl)
     disp(sprintf('%s loaded by GRACE2PLMT',fnpl))
else
    if ~exist(ddir1,'dir')
       error ('The data you asked for are not currently stored.')
       %If you received this error it means the directory datanames does
       %not exist which means you have not downloaded the data or the data
       %was not downloaded to where the program expects 
    end

% DATA CENTER
if strcmp(Pcenter,'GFZ')
    if strcmp(Rlevel,'RL04')
        % Find the coefficient files
        datanames=ls2cell(fullfile(ddir1,'GSM*G---_0004'));
        % Find the error files
        errornames=ls2cell(fullfile(ddir1,'GSM*G---_0004.txt'));
        % Know a priori what the bandwidth of the coefficients is
        Ldata=120;
    elseif strcmp(Rlevel,'RL05')
        % Find the coefficient files
        datanames=ls2cell(fullfile(ddir1,'GSM*G---_005a'));
        % Find the error files
        errornames=ls2cell(fullfile(ddir1,'GSM*G---_005a.txt'));
        % Know a priori what the bandwidth of the coefficients is
        Ldata=90;
    end      
elseif  strcmp(Pcenter,'CSR')
    if strcmp(Rlevel,'RL04')
       datanames=ls2cell(fullfile(ddir1,'GSM*0060_0004'));
       errornames=ls2cell(fullfile(ddir1,'GSM*0060_0004.txt'));
    elseif strcmp(Rlevel,'RL05')
       datanames=ls2cell(fullfile(ddir1,'GSM*0060_0005'));
       % CSR Release Level 5 has no calibrated error files
       %errornames=ls2cell(fullfile(ddir1,'GSM*0060_0005.txt'));
    elseif strcmp(Rlevel,'RL06')
        datanames=ls2cell(fullfile(ddir1,'GSM*BA01_0600'));
        %naming convention was changed for RL06 where BA stands
        %for degree 60 gravity solution and 01 represents unconstrained
        %spherical harmonic solution with a boxcar windowing function 
        %(the meaning of 01 was derived from the L-2 UserHandbook_v4.0).
        %In addition, the BB stands for degree 96 gravity solution with a
        %boxcar windowing function. 
        %The other change was made to the last naming entry, which is now
        %in the form rrvv. In this case rr represents the release number
        %maximum 2 digits and vv represents the maximum 2 digit version 
        %number. So for RL06 the nomenclature is 0600 instead of 0006 for
        %Rl05 previosly. The PID naming convention stays the same.    
    end
   % Know a priori what the bandwidth of the coefficients is
   Ldata=60;
elseif  strcmp(Pcenter,'JPL')
    %Do we want to continue JPL support?
    if strcmp(Rlevel,'RL05')
       datanames=ls2cell(fullfile(ddir1,'GSM*JPLEM*0005'));
       % JPL Release Level 5 has no calibrated error files
       %errornames=ls2cell(fullfile(ddir1,'GSM*0060_0004.txt'));
    else
        error('JPL RL04 solutions not currently stored');
    %elseif strcmp(Rlevel,'RL05');
    %    datanames=ls2cell(fullfile(ddir1,'GSM*JPLEM*005'));
    end
   % Know a priori what the bandwidth of the coefficients is
   Ldata=90;
end

% WGS84 reference SETUP
% For now just hardcode the even zonal coefficients (J), later use
% Frederik's GRS.m program, don't bother with the higher degrees
j2= 0.108262982131e-2*-1.0/(2*2+1)^0.5; % will be row 4
j4=-0.237091120053e-5*-1.0/(2*4+1)^0.5; % will be row 11
% Also useful
a=fralmanac('a_EGM96','Earth');

% C20 CORRECTION SETUP

% Load the C(2,0) coefficients from satellite laser ranging, depending on
% our release level. Note, here the 'NH' part describes a no header version
% of the SLR datafiles.
if strcmp(Rlevel,'RL04')
    slrc20=load(fullfile(getenv('IFILES'),'SLR','C20_RL04_NH.txt'));
elseif strcmp(Rlevel,'RL05')
    slrc20=load(fullfile(getenv('IFILES'),'SLR','C20_RL05_NH.txt'));
elseif strcmp(Rlevel,'RL06')
    slrc20=load(fullfile(getenv('IFILES'),'SLR','C20_RL06_NH.txt'));
   
end
% The sigma error is column 4
slrc20_error=slrc20(:,4)*1e-10;
% Remove the AOD1B model which was removed from the GRACE GSM data but
% restored to this SLR data.  Use the raw value (column 2). See data headers.
slrc20=[slrc20(:,1) slrc20(:,2)-slrc20(:,5)*1e-10];
% Convert the dates to Matlab format
[n,m]=size(slrc20);
slrc20(:,1)=datenum([slrc20(:,1) ones(n,1) ones(n,1)]);
% Make slrc20 relative to the WGS84 ellipsoid
slrc20(:,2) = slrc20(:,2) - j2;

% Degree 1 Correction Setup
if Rlevel=='RL04'
    deg1=load(fullfile(getenv('IFILES'),'GRACE','deg1_RL04_NH.txt'));
elseif Rlevel=='RL05'
    deg1=load(fullfile(getenv('IFILES'),'GRACE','deg1_RL05_NH.txt'));
elseif Rlevel=='RL06'
    deg1=load(fullfile(getenv('IFILES'),'GRACE','deg1_RL06_NH.txt'));
     
end
[n,m] = size(deg1);
dates_str = num2str(deg1(:,1));
deg1dates = datenum([str2num(dates_str(:,1:4)) str2num(dates_str(:,5:6)) 15*ones(n,1)]);
[b,m] = unique(deg1dates);
deg1dates = deg1dates(m);
for i=1:n/2; temp = [deg1(2*i-1,2:7); deg1(2*i,2:7)]; mydeg1(i,:,:) = temp; end;


% Initialize
nmonths = length(datanames);
thedates = zeros(1,nmonths);
[dems,dels]=addmon(Ldata);

% Calibrated errors are normally used instead, but they are kept here for
% completeness.

% Last two columns here are "formal" errors
% l m cosine sine cosine_stddev sine_stddev
potcoffs=nan(nmonths,addmup(Ldata),6);
% Last two columns here are "calibrated" errors
% l m cosine sine
cal_errors=nan(nmonths,addmup(Ldata),4);

% Loop over the months
for index = 1:nmonths 
    % load geopotential coefficients
    fname1=fullfile(ddir1,datanames{index});

    % Open and scan the file (data from all three centers is 10 columns)
    fid = fopen(fname1);
    C = textscan(fid,'%s%s%s%s%s%s%s%s%s%s');
    fclose(fid);

    % Only grab the lines for GRCOF2
    Carray = cat(3,C{:});
    I = strmatch('GRCOF2',Carray(:,1,1),'exact');
    Carray = squeeze(Carray(I,1,:));
    
    % Only want columns 2-7, and as format double
    Carray = Carray(:,2:7);
    lmcosi_month=cellfun(@str2num,Carray);
    % This should be addmup(Ldata)
    [m,n] = size(lmcosi_month);
    
    if strcmp(Pcenter,'GFZ') || strcmp(Pcenter,'CSR')
        % Change the order of the coefficients so that 
        % order m goes as [0 01 012 0123 ...]
        new_ordering = zeros(m,6);
        revdel=[0 Ldata:-1:0];
        i=1;
        for j=1:length(dems)
            k = dels(i)+1 + sum( revdel( (1:dems(i) + 1 ) ) );
            new_ordering(j,:) = lmcosi_month(k,:);    
            i=i+1;
        end
        lmcosi_month = new_ordering;
    elseif strcmp(Pcenter,'JPL')
        % The JPL coefficients are in the order we want already; just need
        % to add the zero and one coefficients
        lmcosi_month = [0 0 0 0 0 0; 1 0 0 0 0 0; 1 1 0 0 0 0; lmcosi_month];
        % Some JPL months are only degree 60, so add in zeros to get up to
        % degree 90
        ldim=length(lmcosi_month);
        if max(lmcosi_month(:,1))==60
            [~,~,~,lmcosiE]=addmon(Ldata);
            lmcosi_month = [lmcosi_month;...
                            lmcosiE(ldim+1:end,1:2) zeros(length(lmcosiE)-ldim,4)];
        end
        
    else
        error('Data center error')
    end
    
    % Remove the mean value of the potential i.e. set 0,0 coff = 0
    lmcosi_month(1,3) = 0;
    
    % Make the geopotential relative to the WGS 84 ellipsoid 
    % A bit redundant since we replace 2,0 shortly
    lmcosi_month(4,3) = lmcosi_month(4,3) - j2;
    lmcosi_month(11,3) = lmcosi_month(11,3) - j4;
    
    % Calculate the midpoint of this data span
    monthstart = datenum([str2num(datanames{index}(7:10))...
                        1 str2num(datanames{index}(11:13))]);
    monthend = datenum([str2num(datanames{index}(15:18))...
                        1 str2num(datanames{index}(19:21))]);
    monthmid = (monthstart+monthend)/2;
    thedates(index) = monthmid;
    
    % Now replace the (2,0) coefficient with the SLR value (referenced to
    % WGS84 above).
    % NOTE: This gives a value different than if you used
    % (column3 - column5) from the SLR data file because that data is
    % referenced to an overall mean, not to the WGS 84 ellipsoid.
    where=slrc20(:,1)>monthstart & slrc20(:,1)<monthend;
    if ~any(where)
        % If there is no SLR value within our specific interval, 
        % use the closest value we have
        [~,where]=min(abs(monthmid - slrc20(:,1)));
    end
    % Need to use slrc20(where,2)
    disp(sprintf('C20 was %12.8e now %12.8e',lmcosi_month(4,3),slrc20(where,2)))
    lmcosi_month(4,3)=slrc20(where,2);
    
    % Now replace the degree 1 coefficients with those from Swenson et al.(2008)
    
    where1=deg1dates(:)>monthstart & deg1dates(:)<monthend;
    if ~any(where1)
        % If there is no Deg1 value within our specific interval, 
        % don't change anything, because we know the first few months are missing
        disp('No change to degree 1')
    else
        disp(['Deg1 value for ' datestr(deg1dates(where1)) ' used.']);
        lmcosi_month(2:3,1:4)=squeeze(mydeg1(where1,:,1:4));
    end
    
    
    % Convert the geopotential coefficients into surface mass density, 
    % if so desired
    if strcmp(units,'SD')
        % Need to make geoid first
        lmcosi_extra=plm2pot([lmcosi_month(:,1:2) lmcosi_month(:,5:6)*a],[],[],[],4);
        lmcosi_month=plm2pot([lmcosi_month(:,1:2) lmcosi_month(:,3:4)*a],[],[],[],4);
        % Add the fornal errors back in columns 5,6
        lmcosi_month=[lmcosi_month lmcosi_extra(:,3:4)];
    end
    
    % Combine into one matrix 
    potcoffs(index,:,:) = lmcosi_month;
    
    %%%
    % CALIBRATED ERRORS
    %%%
    % We have no calibrated errors for CSR release 05, so we have to bypass
    % this section in this case.
if (strcmp(Pcenter,'CSR') || strcmp(Pcenter,'JPL')) && (strcmp(Rlevel,'RL05')...
    || strcmp(Rlevel,'RL06'))
       cal_errors(index,:,:) = [lmcosi_month(:,1:2) zeros(size(lmcosi_month(:,1:2)))];
    else
    fname2=fullfile(ddir1,errornames{index});
    
    % Open and scan the file (data from both Pcenters is 5 columns)
    fid = fopen(fname2);
    E = textscan(fid,'%s%s%s%s%s');
    fclose(fid);

    % Only grab the lines for CALSDV
    Earray = cat(3,E{:});
    I = strmatch('CALSDV',Earray(:,1,1),'exact');
    Earray = squeeze(Earray(I,1,:));
    
    % Only want columns 2-5, and as format double
    Earray = Earray(:,2:5);
    cal_errors_month=cellfun(@str2num,Earray);
    [m,n] = size(cal_errors_month);
    
    % Change the order of the coefficients so that 
    % order m goes as [0 01 012 0123 ...]
    revdel=[0 Ldata:-1:0];
    i=1;
    if strcmp(Pcenter,'CSR')
      new_ordering = zeros(m,4);
      demm=dems;
      for j=1:length(dems)
        k = dels(i)+1 + sum( revdel( (1:dems(i) + 1 ) ) );
        new_ordering(j,:) = cal_errors_month(k,:);
        demm(j)=cal_errors_month(k,2);
        i=i+1;
      end
      cal_errors_month = new_ordering;
    elseif strcmp('GSM-2_2006121-2006151_0028_EIGEN_G---_0004.txt',...
                    errornames(index)) || strcmp(Rlevel,'RL05')
      % for one very odd GFZ file
      new_ordering = zeros(m,4);
      i=4;
      for j=4:length(dems)
        k = dels(i)+1 + sum( revdel( (1:dems(i) + 1 ) ) );
        % This file has only 1 space for the 2,1 coefficients
        if dems(i)==0
          k=k-2;
        else
          k=k-3;
        end
        new_ordering(j-3,:) = cal_errors_month(k,:);
        i=i+1;
      end
      cal_errors_month = new_ordering;
      
    else % for the rest of GFZ, which has slightly less odd formatting
      new_ordering = zeros(m-1,4);
      i=4;
      for j=4:length(dems)
        k = dels(i)+1 + sum( revdel( (1:dems(i) + 1 ) ) );
        % These files have two spaces for the 2,1 coefficients
        if j == 5
          k=k-3;
        else
          k=k-2;
        end
        new_ordering(j-3,:) = cal_errors_month(k,:);
        i=i+1;
      end
      cal_errors_month = new_ordering;
    end
    
    % If from the GFZ data center, add terms for l=0 and 1
    if Pcenter == 'GFZ'
      cal_errors_month = [0 0 0 0; 1 0 0 0; 1 1 0 0; cal_errors_month];
    end
    
    % Replace the C20 error from GRACE with the C20 error from SLR since we
    % used the C20 coefficient from SLR
    disp(sprintf('C20 error was %12.8e now %12.8e',cal_errors_month(4,3),slrc20_error(where)))
    cal_errors_month(4,3)=slrc20_error(where);
    
    % Replace the Deg1 error from GRACE with the Deg1 error from Swenson et al.
    if ~any(where1)
        % Do nothing here
    else
        cal_errors_month(2:3,3:4)=squeeze(mydeg1(where1,:,5:6));
    end
    
    % Convert the geopotential error coefficients into surface mass 
    % density, if so desired
    if strcmp(units,'SD')
        % Need to make geoid first
        a=fralmanac('a_EGM96','Earth');
        cal_errors_month=plm2pot([cal_errors_month(:,1:2) cal_errors_month(:,3:4)*a],[],[],[],4);
    end
    
    % Combine into one matrix
    cal_errors(index,:,:) = cal_errors_month;
    
    end % We have no errors?
    
    disp(['Processed month number: ' num2str(index)]);
end

% Save
save(fnpl,'potcoffs','cal_errors','thedates');

   

end % End if we have a save file already

% Collect output
varns={potcoffs,cal_errors,thedates};
varargout=varns(1:nargout);