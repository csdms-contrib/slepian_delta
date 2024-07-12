% [thedates,GIAt,GIAtU,GIAtL,trend]=CORRECT4GIA(thedates,model,TH,L)
%
% This function accepts an array of dates (in matlab datenum format) and
% calculates the surface mass density change from a certain GIA model at
% each date with reference to the first date.  That is, the GIA models are
% saved as yearly rates of surface mass density change, so there will be a
% different field for each date.
%
% If you want lmcosi fields for GIA, then give only 2 arguments.  If you
% want the GIA projected into a Slepian basis, give that basis.
%
%
% INPUT:
%
% thedates     An array of dates in matlab DATENUM format.  The first date
%                is used as the reference date. [default: monthly points
%                during 2004]
% model        Which GIA model you want.  See references. Options include:
%                'Paulson07'    A model based on the ICE-5G ice load model
%                               of Peltier (2004).  Suitable for both
%                               Antarctica and Greenland.  As corrected by
%                               Geruo A and J. Wahr.  Global.
%                'Wangetal08'   A model based on the older ICE-4G ice
%                               model, and viscosity which varies
%                               laterally.  Suitable for Greenland.
%                'IJ05_R2'      A model based on the Ivins et al. (2013)
%                               ice model.  Suitable for Antarctica.
%                'IJ05'         A model based on the Ivins and James (2005)
%                               ice model.  Suitable for Antarctica.
%                'W12a_v1'      A "best" model from Whitehouse et al (2012)
%                               Suitable only for Antarctica.
%
% TH         Optional [default nothing]:Radius of the concentration
%              region (degrees) OR
%              'england', 'eurasia',  'namerica', 'australia', 'greenland'
%              'africa', 'samerica', 'amazon', 'orinoco', in which case
%              you must have phi,theta,omega all equal to zero OR
%              [lon lat] an ordered list defining a closed curve [degrees]
%              OR a cell containing a region and a buffer such
%              as {'greenland' 0.5}
% Lwindow    Optional [default nothing]: Bandwidth of the window [default
%              if you give a region: bandwidth of the data]
%          If you gave a concentration radius for TH, then you need to
%          specify these:
% phi        Longitude of the center of the tapers (degrees)
% theta      Colatitude of the center of the tapers (degrees)
% omega      Anticlockwise azimuthal rotation of the tapers (degrees)
%
%
%
% OUTPUT:
%
% thedates      Your dates back to you.
% GIAt          If you want SH coefficients for the GIA, then this will be
%                a 3D matrix of GIA fields for each date requested.  The
%                first dimension are the dates.  The second and third
%                dimensions are the familiar lmcosi dimensions
%                (i.e. (L+1)^2 by 4)  Units are kg/m^2
%                OR if you asked for a Slepian projection of this data,
%                then it will be a matrix like "slept" where the first
%                dimension is time and the second dimension is
%                Slepian coefficient.
% GIAtU         Same as GIAt, but if the model has an upper bound
% GIAtL         Same as GIAt, but if the model has a lower bound
% trend         This is the magnitude of the correction in units of Gt/yr.
%                This output will only work if you gave a Slepian basis.
%
%
% NOTES:  It is left to the user to obtain the aforementioned GIA models
% from the authors and save them as appropriate mat files (this function
% assumes rate of change of surface mass density per year).  Perhaps in
% the future these models will be distributed collectively.
%
% SEE ALSO:  PLM2POT
% REFERENCES:
%   Paulson, A., S. Zhong, and J. Wahr. Inference of mantle viscosity from
%    GRACE and relative sea level data, Geophys. J. Int. (2007) 171,
%    497–508. doi: 10.1111/j.1365-246X.2007.03556.x
%
%   Geruo, A., Wahr, J. & Zhong, S. Computations of the viscoelastic response of a
%    3-D compressible Earth to surface loading: An application to Glacial Isostatic
%    Adjustment in Antarctica and Canada. Geophys. J. Int. 192, 557572 (2013).
%
%   Ivins, E. R., T. S. James, J. Wahr, E. J. O. Schrama, F. W. Landerer,
%   and K. M. Simon. Antarctic contribution to sea level rise observed by
%   GRACE with improved GIA correction, Journal of Geophysical Research:
%   Solid Earth, vol. 118, 3126-3141, doi: 10.1002/jgrb.50208, 2013.
%
%   Ivins, E. R., T. S. James. Antarctic glacial isostatic adjustment: a
%   new assessment, Antarctic Science, vol. 17(4), 541-553,
%   doi: 10.1017/S0954102005002968, 2005.
%
%   Wang, H., P. Wu, and W. van der Wal. Using postglacial sea level, crustal
%    velocities and gravity-rate-of-change to constrain the influence of
%    thermal effects on mantle lateral heterogeneities, Journal of
%    Geodynamics (2008) 46, 104-117. doi: 10.1016/j.jog.2008.03.003
%
%   Whitehouse, P.L., Bentley, M.J., Milne, G.A., King, M.A.,
%    Thomas, I.D., 2012. A new glacial isostatic adjustment model for
%    Antarctica: calibrated and tested using observations of relative
%    sea-level change and present-day uplift rates. Geophysical Journal
%    International 190, 1464-1482. doi:10.1111/j.1365-246X.2012.05557.x
%
% Last modified by 
%   charig-at-email.arizona.edu on 5/23/2017
%   williameclee-at-arizona.edu on 7/12/2024

function varargout = correct4gia_new(varargin)
    %% Initialisation
    thedatesDefault = datenum(2004, 1:12, 1); %#ok<DATNM>
    modelDefault = 'W12a_v1';
    DomainDefault = {'greenland', 0.5};
    LDefault = 60;
    phiDefault = 0;
    thetaDefault = 0;
    omegaDefault = 0;
    p = inputParser;
    addOptional(p, 'thedates', thedatesDefault, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'Model', modelDefault, ...
        @(x) ischar(x) || isempty(x));
    addOptional(p, 'Region', DomainDefault, ...
        @(x) isnumeric(x) || iscell(x) || ischar(x));
    addOptional(p, 'L', LDefault, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'phi', phiDefault, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'theta', thetaDefault, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'omega', omegaDefault, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'MoreRegionSpecs', {}, @iscell);
    addParameter(p, 'BeQuiet', false, @islogical);
    parse(p, varargin{:});
    thedates = conddefval(p.Results.thedates, thedatesDefault);
    model = conddefval(p.Results.Model, modelDefault);
    TH = conddefval(p.Results.Region, DomainDefault);
    L = conddefval(p.Results.L, LDefault);
    phi = conddefval(p.Results.phi, phiDefault);
    theta = conddefval(p.Results.theta, thetaDefault);
    omega = conddefval(p.Results.omega, omegaDefault);
    moreRegionSpecs = p.Results.MoreRegionSpecs;
    beQuiet = p.Results.BeQuiet;

    if iscell(TH) && length(TH) > 2
        moreRegionSpecs = {TH{3:end}, moreRegionSpecs{:}}; %#ok<CCAT>
        TH = TH{1:2};
    end

    xver = 0;

    % Where the model save files are kept
    if strncmp(model, 'Morrow', 6)
        defval('ddir1', fullfile(getenv('IFILES'), 'GIA', model(1:6)));
    else
        defval('ddir1', fullfile(getenv('IFILES'), 'GIA', model));
    end

    % And the appropriate name
    fnpl = sprintf('%s/%s_SD.mat', ddir1, model);

    % Load this data (saved as lmcosiM)
    load(fnpl, 'lmcosiM', 'lmcosiU', 'lmcosiL');

    %% Main
    % Convert the model from change per year to change per day
    lmcosiM = [lmcosiM(:, 1:2), lmcosiM(:, 3:4) / 365];

    if exist('lmcosiU', 'var') && exist('lmcosiL', 'var')
        lmcosiU = [lmcosiU(:, 1:2), lmcosiU(:, 3:4) / 365];
        lmcosiL = [lmcosiL(:, 1:2), lmcosiL(:, 3:4) / 365];
    end

    % Reference the date string to the first date
    newdates = thedates - thedates(1);

    % If we have a Slepian basis, do that, if not just do plms
    if exist('TH', 'var')

        if ~exist('L', 'var')
            error('Please specify a basis bandwidth');
        end

        % Figure out if it's lowpass or bandpass
        % lp = length(L) == 1;
        % bp = length(L) == 2;
        maxL = max(L);
        % The spherical harmonic dimension
        % ldim = (L(2 - lp) + 1) ^ 2 - bp * L(1) ^ 2;

        % Project the model
        [~, ~, ~, lmcosiW, ~, ~, ~, ~, ~, ronm] = addmon(maxL);

        if isnumeric(TH)
            [falpha, ~, N, ~, G] = plm2slep_new(lmcosiM, TH, L, phi, theta, omega, "BeQuiet", beQuiet);
        else
            [falpha, ~, N, ~, G] = plm2slep_new(lmcosiM, TH, L, "MoreRegionSpecs", moreRegionSpecs, "BeQuiet", beQuiet);
        end

        if exist('lmcosiU', 'var') && exist('lmcosiL', 'var')
            [falphaU] = plm2slep_new(lmcosiU, TH, L);
            [falphaL] = plm2slep_new(lmcosiL, TH, L);
        end

        % Scale to the new dates
        for i = 1:length(newdates)
            GIAt(i, :) = falpha * newdates(i);

            if exist('lmcosiU', 'var') && exist('lmcosiL', 'var')
                GIAtU(i, :) = falphaU * newdates(i);
                GIAtL(i, :) = falphaL * newdates(i);
            end

            if xver && i == 12
                % Collect the eigenvector output into a format that PLM2XYZ knows how to interpret
                for j = 1:round(N) %#ok<UNRCH>
                    % Create the blanks
                    cosi = lmcosiW(:, 3:4);
                    % Stick in the coefficients of the 1st eigentaper
                    cosi(ronm) = G(:, j);
                    % Construct the full matrix
                    if j == 1
                        CC = [lmcosiW(:, 1:2) cosi * falpha(j) * 365/10];
                    else
                        CC(:, 3:4) = CC(:, 3:4) + cosi * falpha(j) * 365/10;
                    end

                end

                plotplm(CC, [], [], 4, 1)
            end

        end

        % How large is this signal we just made?
        % To answer this we integrate the basis functions and multiply by the
        % coefficients representing the annual rate.
        % Note: the falpha from the models should already be in units of
        % surface mass density change per year.
        for j = 1:round(N)
            cosi = lmcosiW(:, 3:4);
            cosi(ronm) = G(:, j);
            CC{j} = [lmcosiW(:, 1:2), cosi];
        end

        [eigfunINT] = integratebasis(CC, TH, round(N), "MoreRegionSpecs", moreRegionSpecs);
        % Since Int should have units of (fn * m^2), need to go from fractional
        % sphere area to real area.  If the fn is surface density, this output is
        % in kilograms.  Then change the units from kg to Gt in METRIC tons
        eigfunINT = eigfunINT * 4 * pi * 6370000 ^ 2/10 ^ 3/10 ^ 9;
        % functionintegrals = eigfunINT;

        % Now multiply by the appropriate slepcoffs to get the months
        % This becomes alpha by months
        %functimeseries=repmat(eigfunINT',1,nmonths).*sleptdelta(:,1:N)';
        %functimeseries = sleptdelta(:,1:N)';

        total = eigfunINT * (falpha(1:round(N)) * 365); % Back to per year

    else % Just do the plms

        for i = 1:length(newdates)
            GIAt(i, :, :) = [lmcosiM(:, 1:2) lmcosiM(:, 3:4) * newdates(i)];

            if exist('lmcosiU', 'var') && exist('lmcosiL', 'var')
                GIAtU(i, :, :) = [lmcosiU(:, 1:2) lmcosiU(:, 3:4) * newdates(i)];
                GIAtL(i, :, :) = [lmcosiL(:, 1:2) lmcosiL(:, 3:4) * newdates(i)];
            end

            if xver && i == 2
                plotplm(squeeze(lmcosiGIA(i, :, :)), [], [], 4, 1) %#ok<USENS,UNRCH>
                pause
            end

        end

        total = 0;
    end % end if exist

    %% Returning requested outputs
    if exist('lmcosiU', 'var') && exist('lmcosiL', 'var')
        varargout = {thedates, GIAt, GIAtU, GIAtL};
    else
        varargout = {thedates, GIAt, [], [], total};
    end

end
