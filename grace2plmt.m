% [potcoffs,calErrors,thedates]=GRACE2PLMT(Pcenter,Rlevel,units,forcenew)
%
% This program reads in the Level-2 GRACE geoid products from either the CSR or
% GFZ data centers, does some processing, and saves them as a plmt matrix
% in a .mat file.  In particular, the coefficients are reordered to our
% prefered lmcosi format, they are referenced to the WGS84 ellipsoid,
% the C2,0 coefficients are replaced with more accurate measurements from
% satellite laser ranging from Loomis et al, (2020), and the degree one coefficients are
% substituted with those from Sun et al., (2016).  You have the option
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
% Ldata       The bandwidth of the dataproduct that you want [default: 60].
%              In the case where there are more than one product from a
%              datacenter (such as BA 60 or BB 96 standard L2 products)
%              this allows you to choose between them.
% units       'POT' or 'SD' for whether you want geopotential or surface
%               mass density
% forcenew    Whether or not you want to force new generation of a save file
%              (1) or just use the one we already have (0) [default].
%
% OUTPUT:
%
% Returns these variables and saves them in a .mat file:
%    potcoffs       potential coefficients [nmonths x addmup(Ldata) x 4]
%                    these could also be in surface mass density
%    thedates       time stamps in Matlab time
%
% NOTE:
%
%   5/18/2022 Added L to the inputs so that we can utilize more than one
%    product from a data center. A corresponding change has been made in
%    GRACE2SLEPT
%
%   2/19/2021 Formal or calibrated uncertainties have not been reported
%    since RL04 so we will discontinue the output of these. The output
%    variables will be altered so existing scripts will need to be adjusted
%    to take this into account. A corresponding change has been made in
%    GRACE2SLEPT
%
%	GRACE data available from NASA PODAAC:
%	https://podaac.jpl.nasa.gov/
%
%   See gracedeg1.m and gracedeg2.m to see how we update the spherical
%   harmonic degree 1, degree 2, and degree 3 (when available)
%   coefficients. These are usually read from the Technical Notes documents
%   distributed alongside GRACE data.
%
%
%
% EXAMPLE: to make a new save file when you have added more months
% [potcoffs,thedates]=grace2plmt('CSR','RL06',60,'SD',1);
%
% Last modified by charig-at-email.arizona.edu, 5/18/2022
% Last modified by lashokkumar-at-arizona.edu, 11/09/2020
% Last modified by mlubeck-at-email.arizona.edu, 03/18/2019
% Last modified by fjsimons-at-alum.mit.edu, 05/17/2011

% Known problem: the programme will not recalculate the Slepian coefficients when more months of data is added. (WECL, 2024-05-30)

% Determine parameters and set defaults
% Rewritten with inputParser by En-Chi Lee
function varargout = grace2plmt(Pcenter, Rlevel, varargin)
    p = inputParser;
    addRequired(p, 'Pcenter', ...
        @(x) ischar(validatestring(x, {'CSR', 'GFZ', 'JPL'})));
    addRequired(p, 'Rlevel', ...
        @(x) ischar(validatestring(x, {'RL04', 'RL05', 'RL06'})));
    addOptional(p, 'Ldata', 60, ...
        @(x) isnumeric(x) && x > 0);
    addOptional(p, 'units', 'SD', ...
        @(x) ischar(validatestring(x, {'POT', 'SD'})));
    addOptional(p, 'forcenew', 0, ...
        @(x) isnumeric(x) && x == 0 || x == 1);
    parse(p, Pcenter, Rlevel, varargin{:});
    Pcenter = p.Results.Pcenter;
    Rlevel = p.Results.Rlevel;
    Ldata = p.Results.Ldata;
    unit = p.Results.units;
    forcenew = p.Results.forcenew;

    % Where the original data files are kept
    % changed by Will En-Chi Lee for his own directory structure and not to use defval
    if ~isempty(getenv('ORIGINALGRACEDATA'))
        rawDataFolder = fullfile(getenv('ORIGINALGRACEDATA'), ...
            Rlevel, Pcenter);
    elseif ~isempty(getenv('GRACEDATA'))
        rawDataFolder = fullfile(getenv('GRACEDATA'), ...
            'raw', Rlevel, Pcenter);
    else
        rawDataFolder = fullfile(getenv('IFILES'), ...
            'GRACE', 'raw', Rlevel, Pcenter);
    end

    % Where you would like to save the new .mat file
    if ~isempty(getenv('GRACEDATA'))
        processedDataFolder = fullfile(getenv('GRACEDATA'));
    else
        processedDataFolder = fullfile(getenv('IFILES'), 'GRACE');
    end

    % And the name of that save file
    % Rewritten by Will En-Chi Lee to use switch instead of if-else
    switch unit % no otherwise case since input validity is already checked
        case 'SD'
            processedDataFile = sprintf('%s/%s_%s_alldata_%s_%s.mat', ...
                processedDataFolder, Pcenter, Rlevel, num2str(Ldata), unit);
        case 'POT'
            processedDataFile = sprintf('%s/%s_%s_alldata_%s.mat', ...
                processedDataFolder, Pcenter, Rlevel, num2str(Ldata));
    end

    % If this file already exists, load it.  Otherwise, or if we force it, make
    % a new one (e.g. you added extra months to the database).
    if exist(processedDataFile, 'file') == 2 && forcenew == 0
        load(processedDataFile, 'potcoffs', 'thedates')
        fprintf('%s loading %s\n', upper(mfilename), processedDataFile)

        varns = {potcoffs, thedates, thedates};
        varargout = varns(1:nargout);
        return
    end

    if ~exist(rawDataFolder, 'dir') % raw data folder not found
        error ('The data you asked for are not currently stored.')
    end

    % DATA CENTER
    switch Pcenter
        case 'GFZ'

            switch Rlevel
                case 'RL04'
                    % Find the coefficient files
                    dataFiles = ls2cell(fullfile(rawDataFolder, ...
                    'GSM*G---_0004'));
                    % Find the error files
                    errorFiles = ls2cell(fullfile(rawDataFolder, ...
                    'GSM*G---_0004.txt'));
                    % Know a priori what the bandwidth of the coefficients is
                    Ldata = 120;
                case 'RL05'
                    % Find the coefficient files
                    dataFiles = ls2cell(fullfile(rawDataFolder, ...
                    'GSM*G---_005a'));
                    % Find the error files
                    errorFiles = ls2cell(fullfile(rawDataFolder, ...
                    'GSM*G---_005a.txt'));
                    % Know a priori what the bandwidth of the coefficients is
                    Ldata = 90;
            end

        case 'CSR'

            switch Rlevel
                case 'RL04'
                    dataFiles = ls2cell(fullfile(rawDataFolder, ...
                    'GSM*0060_0004'));
                    errorFiles = ls2cell(fullfile(rawDataFolder, ...
                    'GSM*0060_0004.txt'));
                case 'RL05'
                    dataFiles = ls2cell(fullfile(rawDataFolder, ...
                    'GSM*0060_0005'));
                case 'RL06'
                    if Ldata == 60
                        disp(fullfile(rawDataFolder, ...
                        'GSM*BA01_060*'))
                        % dataFiles=[ls2cell(fullfile(ddir1,'GSM*BA01_0600')) ls2cell(fullfile(ddir1,'GSM*BA01_0602'))];
                        dataFiles = ls2cell(fullfile(rawDataFolder, ...
                        'GSM*BA01_060*'));
                    elseif Ldata == 96
                        % dataFiles=[ls2cell(fullfile(ddir1,'GSM*BB01_0600')) ls2cell(fullfile(ddir1,'GSM*BB01_0602'))];
                        dataFiles = ls2cell(fullfile(rawDataFolder, ...
                        'GSM*BB01_060*'));
                    else
                        error(['Solutions with requested L=' num2str(Ldata) ' not currently stored']);
                    end

                    % Naming convention was changed for RL06 where BA stands
                    % for degree 60 gravity solution and 01 represents unconstrained
                    % spherical harmonic solution with a boxcar windowing function
                    % (see L-2 UserHandbook_v4.0).
                    % In addition, the BB stands for degree 96 gravity solution with a
                    % boxcar windowing function.
                    % The other change was made to the last naming entry, which is now
                    % in the form rrvv. In this case rr represents the release number
                    % maximum 2 digits and vv represents the maximum 2 digit version
                    % number. So for RL06 the nomenclature is 0600 instead of 0006 for
                    % Rl05 previosly. The PID naming convention stays the same.
            end

            % Know a priori what the bandwidth of the coefficients is
            % Ldata=60;
        case 'JPL'

            switch Rlevel
                case 'RL05'
                    dataFiles = ls2cell(fullfile(rawDataFolder, ...
                    'GSM*JPLEM*0005'));
                    % JPL Release Level 5 has no calibrated error files
                    %errornames=ls2cell(fullfile(ddir1,'GSM*0060_0004.txt'));
                otherwise
                    error('JPL RL04 solutions not currently stored');
                    %elseif strcmp(Rlevel,'RL05');
                    %    datanames=ls2cell(fullfile(ddir1,'GSM*JPLEM*005'));
            end

            % Know a priori what the bandwidth of the coefficients is
            Ldata = 90;
    end

    % WGS84 reference SETUP
    % For now just hardcode the even zonal coefficients (J), later use
    % Frederik's GRS.m program, don't bother with the higher degrees
    j2 = 0.108262982131e-2 * -1.0 / (2 * 2 + 1) ^ 0.5; % will be row 4
    j4 = -0.237091120053e-5 * -1.0 / (2 * 4 + 1) ^ 0.5; % will be row 11
    % Also useful
    a = fralmanac('a_EGM96', 'Earth');

    % C20 and C30 CORRECTION SETUP
    [slrc20] = gracedeg2(Rlevel);

    % These are not referenced to anything, so make the 2,0 coefficient
    % relative to the WGS84 ellipsoid
    slrc20(:, 2) = slrc20(:, 2) - j2;

    % Degree 1 Correction Setup
    [deg1dates, mydeg1] = gracedeg1(Pcenter, Rlevel);

    % Preallocation
    nmonths = length(dataFiles);
    thedates = zeros(1, nmonths);
    [dems, dels] = addmon(Ldata);
    % Calibrated errors are normally used instead, but they are kept here for
    % completeness.
    % Last two columns here are "formal" errors
    potcoffs = nan(nmonths, addmup(Ldata), 6); % l m cos sin cosStd sinStd
    % Last two columns here are "calibrated" errors
    calErrors = nan(nmonths, addmup(Ldata), 4); % l m cos sin

    %% Loop over the months
    for monthId = 1:nmonths
        disp(['Processing month ' num2str(monthId)]);

        % load geopotential coefficients
        dataFile = fullfile(rawDataFolder, dataFiles{monthId});

        % Open and scan the file (data from all three centers is 10 columns)
        fid = fopen(dataFile);
        C = textscan(fid, '%s%s%s%s%s%s%s%s%s%s');
        fclose(fid);

        % Only grab the lines for GRCOF2
        Carray = cat(3, C{:});
        I = strcmp('GRCOF2', Carray(:, 1, 1));
        Carray = squeeze(Carray(I, 1, :));

        % Only want columns 2-7, and as format double
        Carray = Carray(:, 2:7);
        lmcosi_month = cellfun(@str2num, Carray);
        % This should be addmup(Ldata)
        % [m,n] = size(lmcosi_month);
        m = size(lmcosi_month, 1);

        if strcmp(Pcenter, 'GFZ') || strcmp(Pcenter, 'CSR')
            % Change the order of the coefficients so that
            % order m goes as [0 01 012 0123 ...]
            new_ordering = zeros(m, 6);
            revdel = [0 Ldata:-1:0];
            i = 1;

            for j = 1:length(dems)
                k = dels(i) + 1 + sum(revdel((1:dems(i) + 1)));
                new_ordering(j, :) = lmcosi_month(k, :);
                i = i + 1;
            end

            lmcosi_month = new_ordering;
        elseif strcmp(Pcenter, 'JPL')
            % The JPL coefficients are in the order we want already; just need
            % to add the zero and one coefficients
            lmcosi_month = [0 0 0 0 0 0; 1 0 0 0 0 0; 1 1 0 0 0 0; ...
                                lmcosi_month];
            % Some JPL months are only degree 60, so add in zeros to get up to
            % degree 90
            ldim = length(lmcosi_month);

            if max(lmcosi_month(:, 1)) == 60
                [~, ~, ~, lmcosiE] = addmon(Ldata);
                lmcosi_month = [lmcosi_month; ...
                                    lmcosiE(ldim + 1:end, 1:2) zeros(length(lmcosiE) - ldim, 4)];
            end

        end

        % Remove the mean value of the potential i.e. set 0,0 coff = 0
        lmcosi_month(1, 3) = 0;

        % Make the geopotential relative to the WGS 84 ellipsoid
        % A bit redundant since we replace 2,0 shortly
        lmcosi_month(4, 3) = lmcosi_month(4, 3) - j2;
        lmcosi_month(11, 3) = lmcosi_month(11, 3) - j4;

        % Calculate the midpoint of this data span
        monthstart = datenum([str2num(dataFiles{monthId}(7:10)) ...
                                  1 str2num(dataFiles{monthId}(11:13))]);
        monthend = datenum([str2num(dataFiles{monthId}(15:18)) ...
                                1 str2num(dataFiles{monthId}(19:21))]);
        monthmid = (monthstart + monthend) / 2;
        thedates(monthId) = monthmid;

        disp(['Month midpoint date: ' ...
                  char(datetime(thedates(monthId), ...
                  'ConvertFrom', 'datenum'))]);

        % Now replace the (2,0) and (3,0) (when available) coefficient with
        % the SLR value (referenced to WGS84 above).
        %     % NOTE: This gives a value different than if you used
        %     % (column3 - column5) from the SLR data file because that data is
        %     % referenced to an overall mean, not to the WGS 84 ellipsoid.
        %     %if index==111; keyboard; end;
        where = slrc20(:, 1) > monthstart & slrc20(:, 1) < monthend;

        if ~any(where)
            % If there is no SLR value within our specific interval,
            % use the closest value we have
            [~, where] = min(abs(monthmid - slrc20(:, 1)));
            disp('No SLR coefficients for this month, used nearest available.')
        elseif nnz(where) > 1
            % We have more than one month. Use the closest value
            [~, where] = min(abs(monthmid - slrc20(:, 1)));
            disp(['More than one degree 2 value detected, likely from overlapping months. Used nearest: ' datestr(slrc20(where, 1))])
        end

        % else we just use "where" which only has 1 valid entry

        % Need to use slrc20(where,2)
        fprintf('C20: %+10.6e -> %+10.6e\n', ...
            lmcosi_month(4, 3), slrc20(where, 2))
        lmcosi_month(4, 3) = slrc20(where, 2);
        % Now replace C3,0 if possible
        if isfinite(slrc20(where, 3))
            % Here we have a C3,0 SLR value available, so substitute.
            fprintf('C30: %+10.6e -> %+10.6e\n', ...
                lmcosi_month(7, 3), slrc20(where, 3))
            lmcosi_month(7, 3) = slrc20(where, 3);
        end

        % Now replace the degree 1 coefficients with those from GRACEDEG1.m and
        % references therin.

        % Find a degree 1 data point that is within the month of our GRACE data
        where1 = deg1dates(:) > monthstart & deg1dates(:) < monthend;

        if ~any(where1)
            % If there is no Deg1 value within our specific interval,
            % don't change anything, because we know a few months are missing
            disp('No change to degree 1')
        elseif nnz(where1) > 1
            % We have more than one month. Use the closest value
            disp('More than one degree 1 value detected')
            [~, where1] = min(abs(monthmid - deg1dates));
            lmcosi_month(2:3, 1:4) = squeeze(mydeg1(where1, :, 1:4));
            disp(['Deg1 value for ' datestr(deg1dates(where1)) ' used.']);
        else
            disp(['Deg1 value for ' datestr(deg1dates(where1)) ' used.']);
            lmcosi_month(2:3, 1:4) = squeeze(mydeg1(where1, :, 1:4));
        end

        % Convert the geopotential coefficients into surface mass density,
        % if so desired
        if strcmp(unit, 'SD')
            % Need to make geoid first
            lmcosi_extra = plm2pot([lmcosi_month(:, 1:2) ...
                                        lmcosi_month(:, 5:6) * a], [], [], [], 4);
            lmcosi_month = plm2pot([lmcosi_month(:, 1:2) ...
                                        lmcosi_month(:, 3:4) * a], [], [], [], 4);
            % Add the fornal errors back in columns 5,6
            lmcosi_month = [lmcosi_month lmcosi_extra(:, 3:4)];
        end

        % Combine into one matrix
        potcoffs(monthId, :, :) = lmcosi_month;

        %%%
        % CALIBRATED ERRORS
        %%%
        % We have no calibrated errors for CSR release 05, so we have to bypass
        % this section in this case.
        if (strcmp(Pcenter, 'CSR') || strcmp(Pcenter, 'JPL')) && (strcmp(Rlevel, 'RL05') ...
                || strcmp(Rlevel, 'RL06'))
            calErrors(monthId, :, :) = [lmcosi_month(:, 1:2) zeros(size(lmcosi_month(:, 1:2)))];
        else
            fname2 = fullfile(rawDataFolder, errorFiles{monthId});

            % Open and scan the file (data from both Pcenters is 5 columns)
            fid = fopen(fname2);
            E = textscan(fid, '%s%s%s%s%s');
            fclose(fid);

            % Only grab the lines for CALSDV
            Earray = cat(3, E{:});
            I = strcmp('CALSDV', Earray(:, 1, 1));
            Earray = squeeze(Earray(I, 1, :));

            % Only want columns 2-5, and as format double
            Earray = Earray(:, 2:5);
            cal_errors_month = cellfun(@str2num, Earray);
            % [m,n] = size(cal_errors_month);
            m = size(cal_errors_month, 1);

            % Change the order of the coefficients so that
            % order m goes as [0 01 012 0123 ...]
            revdel = [0 Ldata:-1:0];
            i = 1;

            if strcmp(Pcenter, 'CSR')
                new_ordering = zeros(m, 4);
                demm = dems;

                for j = 1:length(dems)
                    k = dels(i) + 1 + sum(revdel((1:dems(i) + 1)));
                    new_ordering(j, :) = cal_errors_month(k, :);
                    demm(j) = cal_errors_month(k, 2);
                    i = i + 1;
                end

                cal_errors_month = new_ordering;
            elseif strcmp('GSM-2_2006121-2006151_0028_EIGEN_G---_0004.txt', ...
                    errorFiles(monthId)) || strcmp(Rlevel, 'RL05')
                % for one very odd GFZ file
                new_ordering = zeros(m, 4);
                i = 4;

                for j = 4:length(dems)
                    k = dels(i) + 1 + sum(revdel((1:dems(i) + 1)));
                    % This file has only 1 space for the 2,1 coefficients
                    if dems(i) == 0
                        k = k - 2;
                    else
                        k = k - 3;
                    end

                    new_ordering(j - 3, :) = cal_errors_month(k, :);
                    i = i + 1;
                end

                cal_errors_month = new_ordering;

            else % for the rest of GFZ, which has slightly less odd formatting
                new_ordering = zeros(m - 1, 4);
                i = 4;

                for j = 4:length(dems)
                    k = dels(i) + 1 + sum(revdel((1:dems(i) + 1)));
                    % These files have two spaces for the 2,1 coefficients
                    if j == 5
                        k = k - 3;
                    else
                        k = k - 2;
                    end

                    new_ordering(j - 3, :) = cal_errors_month(k, :);
                    i = i + 1;
                end

                cal_errors_month = new_ordering;
            end

            % If from the GFZ data center, add terms for l=0 and 1
            % if Pcenter == 'GFZ'
            if strcmp(Pcenter, 'GFZ')
                cal_errors_month = [0 0 0 0; 1 0 0 0; 1 1 0 0; cal_errors_month];
            end

            % Replace the C20 error from GRACE with the C20 error from SLR since we
            % used the C20 coefficient from SLR
            fprintf('C20 error was %12.8e now %12.8e\n', cal_errors_month(4, 3), slrc20_error(where))
            cal_errors_month(4, 3) = slrc20_error(where);

            % Replace the Deg1 error from GRACE with the Deg1 error from Swenson et al.
            if ~any(where1)
                % Do nothing here
            else
                cal_errors_month(2:3, 3:4) = squeeze(mydeg1(where1, :, 5:6));
            end

            % Convert the geopotential error coefficients into surface mass
            % density, if so desired
            if strcmp(unit, 'SD')
                % Need to make geoid first
                a = fralmanac('a_EGM96', 'Earth');
                cal_errors_month = plm2pot([cal_errors_month(:, 1:2) cal_errors_month(:, 3:4) * a], [], [], [], 4);
            end

            % Combine into one matrix
            calErrors(monthId, :, :) = cal_errors_month;

        end % We have no errors?

        disp(' ')
    end

    % Save
    save(processedDataFile, 'potcoffs', 'thedates');

    % Collect output
    % Here we have "thedates" twice so that we don't break older code. But in
    % the future we will fix this so that we only have the two output
    varns = {potcoffs, thedates, thedates};
    varargout = varns(1:nargout);
end
