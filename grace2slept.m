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
function varargout = grace2slept(varargin)
    %% Parsing inputs and initialisation
    dataProductDefault = {'CSR', 'RL06', 60};
    regionDefault = 'greenland';
    bufferDefault = 0;
    LwindowDefault = 18;
    taperDefault = 0;
    JDefault = [];
    unitsDefault = 'POT';
    forcenewDefault = 0;

    p = inputParser;
    addOptional(p, 'DataProduct', dataProductDefault, ...
        @(x) iscell(x) && length(x) == 3 ...
        && ischar(x{1}) && ischar(x{2}) && isnumeric(x{3}));
    addOptional(p, 'TH', regionDefault, ...
        @(x) (ischar(x)) || ...
        (isnumeric(x) && size(x, 2) == 2) || ...
        (isempty(x)));
    addOptional(p, 'XY_buffer', bufferDefault, ...
        @(x) (isnumeric(x) && x >= 0) || (isempty(x)));
    addOptional(p, 'Lwindow', LwindowDefault, ...
        @(x) (isnumeric(x) && isscalar(x) && x > 0) || (isempty(x)));
    addOptional(p, 'phi', taperDefault, ...
        @(x) (isnumeric(x) && isscalar(x)) || (isempty(x)));
    addOptional(p, 'theta', taperDefault, ...
        @(x) (isnumeric(x) && isscalar(x)) || (isempty(x)));
    addOptional(p, 'omega', taperDefault, ...
        @(x) (isnumeric(x) && isscalar(x)) || (isempty(x)));
    addOptional(p, 'J', JDefault, ...
        @(x) ((isnumeric(x) && isscalar(x) && x > 0) || strcmp(x, 'N')) || ...
        (isempty(x)));
    addOptional(p, 'units', unitsDefault, ...
        @(x) (ischar(validatestring(x, {'POT', 'SD'}))) || (isempty(x)));
    addOptional(p, 'forcenew', forcenewDefault, ...
        @(x) (isnumeric(x) && (x == 0 || x == 1)) || (isempty(x)));
    addParameter(p, 'upscale', 0, @(x) isnumeric(x) && x >= 0);
    parse(p, varargin{:});

    dataProduct = conddefval(p.Results.DataProduct, dataProductDefault);
    region = conddefval(p.Results.TH, regionDefault);
    buffer = conddefval(p.Results.XY_buffer, bufferDefault);
    Lwindow = conddefval(p.Results.Lwindow, LwindowDefault);
    phi = conddefval(p.Results.phi, taperDefault);
    theta = conddefval(p.Results.theta, taperDefault);
    omega = conddefval(p.Results.omega, taperDefault);
    unit = conddefval(p.Results.units, unitsDefault);
    forcenew = conddefval(p.Results.forcenew, forcenewDefault);
    upscale = p.Results.upscale;
    % note: J is defined later in the code
    % pars and inout not used
    % pars = 10;
    % inout = 'out';

    dataCentre = dataProduct{1};
    releaseLevel = dataProduct{2};
    Ldata = dataProduct{3};
    dataproductstring = [dataCentre releaseLevel num2str(Ldata)];

    % Figure out if it's lowpass or bandpass
    lp = isscalar(Lwindow);
    bp = length(Lwindow) == 2;
    maxL = max(Lwindow);
    % The spherical harmonic dimension
    ldim = (Lwindow(2 - lp) + 1) ^ 2 - bp * Lwindow(1) ^ 2;
    J = conddefval(p.Results.J, ldim);

    %% Where you would like to save the new .mat file
    if ~isempty(getenv("GRACE"))
        outputFolder = fullfile(getenv("GRACE"), 'SlepianExpansions');
    else
        outputFolder = fullfile(getenv("IFILES"), 'GRACE', 'SlepianExpansions');
    end

    % Get the remaining file names
    if ~isstr(region) && length(region) == 1 % POLAR CAPS
        % Check if you want the Shannon number of eigenfunctions
        if strcmp(J, 'N')
            J = round((Lwindow + 1) ^ 2 * spharea(region));
        end

        if lp
            fnpl = sprintf('%s/grace2slept-%s-CAP-%i-%i-%i-%s.mat', ...
                outputFolder, dataproductstring, region, Lwindow, J, unit);
        elseif bp
            fnpl = sprintf('%s/grace2sleptbl-%s-CAP-%i-%i-%i-%i-%s.mat', ...
                outputFolder, dataproductstring, region, Lwindow(1), Lwindow(2), J, unit);
        else
            error('The degree range is either one or two numbers')
        end

    else % GEOGRAPHICAL REGIONS and XY REGIONS
        % Check if you want the Shannon number of eigenfunctions
        if strcmp(J, 'N')
            J = round((Lwindow + 1) ^ 2 * spharea(region));
        end

        if isstr(region) % Geographic
            % Here, TH gets passed to glmalpha, and glmalpha will interpret
            % either the cell of the region
            if buffer ~= 0
                region = {region buffer};
                h = [region{1} num2str(buffer)];
            else
                h = region;
            end

        else % Closed coordinates (make a hash)
            h = hash(region, 'sha1');
        end

        % The name of the save file
        if lp
            fnpl = sprintf('%s/grace2slept-%s-%s-%i-%i-%s.mat', ...
                outputFolder, dataproductstring, h, Lwindow, J, unit);
        elseif bp
            fnpl = sprintf('%s/grace2sleptbl-%s-%s-%i-%i-%i-%s.mat', ...
                outputFolder, dataproductstring, h, Lwindow(1), Lwindow(2), J, unit);
        else
            error('The degree range is either one or two numbers')
        end

    end

    % GET THE SLEPIAN BASIS WE WANT
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
        [G, V, ~, ~, N, ~, MTAP, ~] ...
            = glmalpha(region, Lwindow, upscale, 0, [], [], J);
    else
        % Need to get a complete GLMALPHA but for the rotated basis
        % Definitely, "single-order" has lost its meaning here, but the MTAP
        % will still identify what the order of the unrotated original was
        [G, V, ~, ~, N, ~, MTAP, ~] ...
            = glmalphapto(region, Lwindow, phi, theta, omega);
        % Since GLMALPHAPTO currently has no option to limit a basis to J, do it here
        G = G(:, 1:J);
    end

    % Sort by decreasing eigenvalue
    [V, vi] = sort(V, 'descend');
    G = G(:, vi);

    if ~isnan(MTAP)
        MTAP = MTAP(vi);
    end

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
    if exist(fnpl, 'file') == 2 && forcenew == 0
        load(fnpl)
        fprintf('%s loading %s\n', upper(mfilename), fnpl)
    else
        % Use GRACE2PLMT to get the GRACE data.  This way, if we have it saved,
        % there is no need to scan the month files again.  GRACE2PLMT takes
        % care of the WGS84 adjustment, the degree 1 correction, and the C20 correction.
        [potcoffs, cal_errors, thedates] = ...
            grace2plmt(dataCentre, releaseLevel, Ldata, unit, 0);
        % *** Here I changed this. Run grace2plmt once to update your data, and
        % then when you call forcenew=1 from now on it will just update the
        % expansion

        % Initialize new coefficients
        nmonths = length(thedates);
        slepcoffs = nan(nmonths, J);
        slepcalerrors = nan(nmonths, J);

        % Limit everything to the window bandwidth
        potcoffsW = potcoffs(:, 1:size(lmcosiW, 1), 1:4);
        %cal_errorsW = cal_errors(:,1:size(lmcosiW,1),1:4);

        % Loop over the months
        for monthId = 1:nmonths
            % Expand this month's POTENTIAL into the Slepian basis
            potcoffs = squeeze(potcoffsW(monthId, :, :));
            slepcoffs(monthId, :) = ...
                potcoffs(2 * size(potcoffs, 1) + ...
                ronmW(1:(maxL + 1) ^ 2))' * G;

            % Expand this month of CALIBRATED ERRORS into the Slepian basis
            %calerrors_month=squeeze(cal_errorsW(index,:,:));
            %slepcalerrors(index,:) = ...
            %    calerrors_month(2*size(calerrors_month,1)+ronmW(1:(maxL+1)^2))'*G;
        end

        % SAVE
        % Here we have "thedates" twice so that we don't break older code.
        % But in the future we will fix this so that we don't have cal
        % errors
        % Don't save the kernel and eigenfunctions because we already have
        % this info saved and can load from GLMALPHA
        % save(fnpl,'slepcoffs','calerrors','thedates','G','CC','V');
        save(fnpl, 'slepcoffs', 'thedates', 'thedates');

    end % End if we have a save file already

    % Collect output
    varns = {slepcoffs, thedates, thedates, region, G, CC, V, N};
    varargout = varns(1:nargout);
