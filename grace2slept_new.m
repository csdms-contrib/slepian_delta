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
%   TH - The region
%       If there was buffering, this will be a XY array of coordinates, which you can use with SPHAREA to get the Shannon number.
%  G            The unitary matrix of localization coefficients
%  CC           A cell array with cosine/sine coefficients eigenfunctions
%  V            The eigenvalues in this ordering
%  N            The Shannon number
%
% SEE ALSO: PLM2SLEP
%
% Last modified by charig-at-princeton.edu, 05/18/2022
% Last modified by fjsimons-at-alum.mit.edu, 06/26/2012

function varargout = grace2slept_new(varargin)
    %% Initialisation
    % Parse inputs
    [~, region, buf, moreRegionSpecs, Lwindow, ...
         phi, theta, omega, unit, upscale, ...
         dataCentre, releaseLevel, Ldata, dataProductString, truncation, ...
         forceReload, saveData, beQuiet] = ...
        parseinputs(varargin);

    % Figure out if it's low-pass or band-pass
    lp = isscalar(Lwindow);
    bp = length(Lwindow) == 2;

    if ~(lp || bp)
        error('The degree range is either one or two numbers')
    end

    maxL = max(Lwindow);
    % The spherical harmonic dimension
    ldim = (Lwindow(2 - lp) + 1) ^ 2 - bp * Lwindow(1) ^ 2;

    % Check if you want the Shannon number of eigenfunctions
    if strcmp(truncation, 'N')
        truncation = round((Lwindow + 1) ^ 2 * spharea(region));
    else
        truncation = conddefval(truncation, ldim);
    end

    % Output file
    [outputFilePath, region] = getoutputfile(region, buf, moreRegionSpecs, Lwindow, ...
        dataProductString, truncation, unit, lp, bp);

    %% Constructing the Slepian basis
    [~, ~, ~, lmcosiW, ~, ~, ~, ~, ~, ronmW] = addmon(maxL);
    % NOTE: a) Could have used PLM2SLEP but since we need the same operation for many
    %         months, slightly better to load the Slepian basis once, and repeatedly
    %         multiply by G.  Mostly copied from PLM2SLEP.
    %       b) The kernel and eigenfunctions are large to save for large Lwindow.
    %          So instead just load them from either GLMALPHA or GLMALPHAPTO.
    %       c) GLMALPHA now handles bandpass

    % If it is the standard North-Polar cap or a geographic region, it's easy
    if phi == 0 && theta == 0 && omega == 0
        % Get the Slepian basis; definitely not block-sorted as for the rotated
        % versions this will make no sense at all anymore
        moreRegionSpecs = {"Upscale", upscale, moreRegionSpecs{:}}; %#ok<CCAT>

        if iscell(region)

            if length(region) >= 3
                [G, V, ~, ~, N] ...
                    = glmalpha_new(region, Lwindow, 'J', truncation, "BeQuiet", beQuiet);
            else
                [G, V, ~, ~, N] ...
                    = glmalpha_new(region{1}, Lwindow, moreRegionSpecs, 0, 'J', truncation, "BeQuiet", beQuiet);
            end

        else
            [G, V, ~, ~, N] ...
                = glmalpha_new(region, Lwindow, moreRegionSpecs, 0, 'J', truncation, "BeQuiet", beQuiet);
        end

    else
        % Need to get a complete GLMALPHA but for the rotated basis
        % Definitely, "single-order" has lost its meaning here, but the MTAP
        % will still identify what the order of the unrotated original was
        [G, V, ~, ~, N] ...
            = glmalphapto(region, Lwindow, phi, theta, omega);
        % Since GLMALPHAPTO currently has no option to limit a basis to J, do it here
        G = G(:, 1:truncation);
    end

    % Sort by decreasing eigenvalue
    [V, vi] = sort(V, 'descend');
    G = G(:, vi);

    % If you don't do this, the eigenfunctions are ordered in the way
    %   that they correspond to single-orders back when, unrotated, they
    %   belonged to a polar cap, and the eigenvalues are sorted within
    %   these blocks. This is useful for, e.g. SPIE2009_1 a la SDSNEEUW.
    % Collect the eigenvector output into a format that PLM2XYZ knows how to interpret
    CC = cell(1, size(G, 2));

    for j = 1:size(G, 2)
        % Create the blanks
        cosi = lmcosiW(:, 3:4);
        % Stick in the coefficients of the 1st eigentaper
        cosi(ronmW) = G(:, j);
        % Construct the full matrix
        CC{j} = [lmcosiW(:, 1:2) cosi];
    end

    % INITILIZATION COMPLETE

    % If this expansion already exists, load it.  Otherwise, or if we force
    % it, make a new one (e.g. if you added extra months to the database).
    if exist(outputFilePath, 'file') == 2 && ~forceReload
        load(outputFilePath, 'slepcoffs', 'thedates')

        if ~beQuiet
            fprintf('%s loading %s\n', upper(mfilename), outputFilePath)
        end

        varargout = {slepcoffs, thedates, thedates, region, G, CC, V, N};

        return
    end

    % Use GRACE2PLMT to get the GRACE data.  This way, if we have it saved,
    % there is no need to scan the month files again.  GRACE2PLMT takes
    % care of the WGS84 adjustment, the degree 1 correction, and the C20 correction.
    [sphCoeffs, cal_errors, thedates] = ...
        grace2plmt(dataCentre, releaseLevel, Ldata, unit, 0); %#ok<ASGLU>
    % *** Here I changed this. Run grace2plmt once to update your data, and
    % then when you call forcenew=1 from now on it will just update the
    % expansion

    % Initialize new coefficients
    nmonths = length(thedates);
    slepcoffs = nan(nmonths, truncation);
    % slepcalerrors = nan(nmonths, truncation);

    % Limit everything to the window bandwidth
    potcoffsW = sphCoeffs(:, 1:size(lmcosiW, 1), 1:4);
    %cal_errorsW = cal_errors(:,1:size(lmcosiW,1),1:4);

    % Loop over the months
    for monthId = 1:nmonths
        % Expand this month's POTENTIAL into the Slepian basis
        sphCoeffs = squeeze(potcoffsW(monthId, :, :));
        slepcoffs(monthId, :) = ...
            sphCoeffs(2 * size(sphCoeffs, 1) + ...
            ronmW(1:(maxL + 1) ^ 2))' * G;

        % Expand this month of CALIBRATED ERRORS into the Slepian basis
        %calerrors_month=squeeze(cal_errorsW(index,:,:));
        %slepcalerrors(index,:) = ...
        %    calerrors_month(2*size(calerrors_month,1)+ronmW(1:(maxL+1)^2))'*G;
    end

    if saveData
        save(outputFilePath, 'slepcoffs', 'thedates', 'thedates');

        if ~beQuiet
            fprintf('%s saving %s\n', upper(mfilename), outputFilePath)
        end

    end

    % Collect output
    varargout = {slepcoffs, thedates, thedates, region, G, CC, V, N};

end

%% Subfunctions
function varargout = parseinputs(vArargin)
    dataProductDefault = {'CSR', 'RL06', 60};
    regionDefault = 'greenland';
    bufDefault = 0;
    LwindowDefault = 18;
    taperDefault = 0;
    JDefault = [];
    unitsDefault = 'POT';
    forceReloadDefault = false;

    p = inputParser;
    addOptional(p, 'DataProduct', dataProductDefault, ...
        @(x) iscell(x) && length(x) == 3 ...
        && ischar(x{1}) && ischar(x{2}) && isnumeric(x{3}));
    addOptional(p, 'Region', regionDefault, ...
        @(x) (ischar(x)) || isstring(x) || iscell(x) || ...
        (isnumeric(x) && size(x, 2) == 2) || ...
        (isempty(x)));
    addOptional(p, 'Buffer', [], ...
        @(x) (isnumeric(x) && x >= 0) || isempty(x));
    addOptional(p, 'Lwindow', LwindowDefault, ...
        @(x) (isnumeric(x) && (length(x) <= 2)) || (isempty(x)));
    addOptional(p, 'phi', taperDefault, ...
        @(x) (isnumeric(x) && isscalar(x)) || (isempty(x)));
    addOptional(p, 'theta', taperDefault, ...
        @(x) (isnumeric(x) && isscalar(x)) || (isempty(x)));
    addOptional(p, 'omega', taperDefault, ...
        @(x) (isnumeric(x) && isscalar(x)) || (isempty(x)));
    addOptional(p, 'Truncation', JDefault, ...
        @(x) ((isnumeric(x) && isscalar(x) && x > 0) || ...
        strcmp(x, 'N')) || (isempty(x)));
    addOptional(p, 'Unit', unitsDefault, ...
        @(x) (ischar(validatestring(x, {'POT', 'SD'}))) || (isempty(x)));
    addOptional(p, 'ForceNew', forceReloadDefault, ...
        @(x) (isnumeric(x) && (x == 0 || x == 1)) || islogical(x) ...
        || (isempty(x)));
    addOptional(p, 'MoreRegionSpecs', {}, @iscell);
    addOptional(p, 'Upscale', 0, @(x) isnumeric(x) && x >= 0);
    addParameter(p, 'BeQuiet', false, @islogical);
    addParameter(p, 'Save', true, @islogical);
    parse(p, vArargin{:});

    dataProduct = conddefval(p.Results.DataProduct, dataProductDefault);
    region = conddefval(p.Results.Region, regionDefault);

    if iscell(region) && length(region) >= 2
        bufDefault = region{2};
    end

    buf = conddefval(p.Results.Buffer, bufDefault);
    Lwindow = conddefval(p.Results.Lwindow, LwindowDefault);
    phi = conddefval(p.Results.phi, taperDefault);
    theta = conddefval(p.Results.theta, taperDefault);
    omega = conddefval(p.Results.omega, taperDefault);
    unit = conddefval(p.Results.Unit, unitsDefault);
    upscale = p.Results.Upscale;
    moreRegionSpecs = p.Results.MoreRegionSpecs;
    J = conddefval(p.Results.Truncation, JDefault);
    forcenew = conddefval(logical(p.Results.ForceNew), forceReloadDefault);
    beQuiet = p.Results.BeQuiet;
    saveData = p.Results.Save;

    dataCentre = dataProduct{1};
    releaseLevel = dataProduct{2};
    Ldata = dataProduct{3};
    dataproductstring = [dataCentre releaseLevel num2str(Ldata)];

    if iscell(region) && length(region) >= 2

        if length(region) > 2
            moreRegionSpecs = {region{3:end}, moreRegionSpecs{:}}; %#ok<CCAT>
        end

        % region = region{1};
    end

    varargout = ...
        {dataProduct, region, buf, moreRegionSpecs, Lwindow, ...
         phi, theta, omega, unit, upscale, ...
         dataCentre, releaseLevel, Ldata, dataproductstring, J, ...
         forcenew, saveData, beQuiet};
end

function [outputFilePath, region] = getoutputfile(region, buf, moreRegionSpecs, Lwindow, ...
        productString, truncation, unit, lp, bp)

    % Folder
    if ~isempty(getenv('GRACE'))
        outputFolder = fullfile(getenv('GRACE'), ...
        'SlepianExpansions');
    else
        outputFolder = fullfile(getenv('IFILES'), ...
            'GRACE', 'SlepianExpansions');
    end

    % File name
    switch regiontype(region)
        case 'polar'

            if lp
                outputFileName = sprintf( ...
                    'grace2slept-%s-CAP-%i-%i-%i-%s.mat', ...
                    productString, region, ...
                    Lwindow, truncation, unit);
            elseif bp
                outputFileName = sprintf( ...
                    'grace2sleptbl-%s-CAP-%i-%i-%i-%i-%s.mat', ...
                    productString, region, ...
                    Lwindow(1), Lwindow(2), truncation, unit);
            end

        case {'geographical', 'XY'}

            switch regiontype(region)
                case 'geographical'
                    % Here, TH gets passed to glmalpha, and glmalpha will interpret
                    % either the cell of the region
                    if iscell(region) && length(region) >= 3
                        regionString = [region{1}, dataattrchar("Buffer", region{3:end})];
                    else
                        regionString = [region, dataattrchar('Buffer', buf, moreRegionSpecs{:})];

                        if buf ~= 0
                            region = {region buf};
                        end

                    end

                    % else % Closed coordinates (make a hash)
                case 'XY'
                    regionString = hash(region, 'sha1');
            end

            % The name of the save file
            if lp
                outputFileName = sprintf( ...
                    'grace2slept-%s-%s-%i-%i-%s.mat', ...
                    productString, regionString, ...
                    Lwindow, truncation, unit);
            elseif bp
                outputFileName = sprintf( ...
                    'grace2sleptbl-%s-%s-%i-%i-%i-%i-%s.mat', ...
                    productString, regionString, ...
                    Lwindow(1), Lwindow(2), truncation, unit);
            end

    end

    outputFilePath = fullfile(outputFolder, outputFileName);

end

function regionType = regiontype(region)

    if ischar(region) || isstring(region) || iscell(region)
        regionType = 'geographical';
    elseif isnumeric(region) && isscalar(region)
        regionType = 'polar';
    elseif isnumeric(region) && ~isscalar(region)
        regionType = 'XY';
    else
        error('Unknown region type')
    end

end
