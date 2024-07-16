% [eigfunINT]=INTEGRATEBASIS(CC,TH,J)
%
% Accepts a Slepian basis and integrates the functions within a region.
% The integrals reported are in fractional sphere area.  You may want to
% convert this to real sphere area via x(4*pi*radius^2).
%
% INPUT:
%
% CC      The eigenfunctions.  This can be either the G you get from
%          GLMALPHA or a cell array of the functions as in LOCALIZATION.
% TH      The region.  Can be formatted multiple ways as in GLMALPHA
% J       How many functions you want to do this for [DEFAULT: all]
% phi     Longitude of the center (degrees)
% theta   Colatitude of the center (degrees)
%           Note: there is no omega here because a rotation of the
%           polar cap functions does not affect their integration
%
% OUTPUT:
%
% eigfunINT    The integrals of the eigenfunctions over the region.  These
%               are in real sphere area, depending on your radius.
%
%
% SEE ALSO: PLM2AVG, PLM2AVGP
%
%
% Last modified by charig-at-email.arizona.edu, 11/01/2016

function varargout = integratebasis_new(varargin)
    [CC, TH, J, phi, theta, MoreRegionSpecs] = parseinputs(varargin);

    %% Initialisation
    % Sort out what CC is
    if isstr(CC)
        % Evaluate the specified expression
        eval(CC);
    end

    if isfloat(CC)
        % We must have a G matrix from glmalpha (already sorted)
        % Reorder them into a cell array
        [n, m] = size(CC);
        defval('J', m);
        L = sqrt(n) - 1;
        % This should be an integer
        if (floor(L) ~= L)
            error('Something fishy about your L');
        end

        [~, ~, ~, lmcosi, ~, ~, ~, ~, ~, ronm] = addmon(L);
        % Collect the eigenvector output into a format that PLM2XYZ knows how to interpret
        for j = 1:J
            % Create the blanks
            cosi = lmcosi(:, 3:4);
            % Stick in the coefficients of the 1st eigentaper
            cosi(ronm) = CC(:, j);
            % Construct the full matrix
            CC2{j} = [lmcosi(:, 1:2) cosi];
        end

        CC = CC2;
    elseif iscell(CC)
        % We are all good, just get the size.
        m = size(CC, 2);
        defval('J', m);
        defval('L', addmup(size(CC{1}, 1), 'r'));
    else
        error('What format is CC?');
    end

    %% Main
    switch regiontype(TH)
        case 'geographic'
            % Default parameters
            defval('buf', 0);
            defval('pars', 0);
            defval('rotb', 0);

            if iscell(TH)
                dom = TH{1};
                buf = TH{2};
            else
                dom = TH;
            end

            % Rotate if needed
            try
                rotb = feval(dom, 'rotated');
            catch
            end

            if length(rotb) == 1 && rotb
                % Return the coordinates and do the rotation
                [XY, lonc, latc] = feval(dom, pars, buf);
                [thetap, phip, ~] = ...
                    rottp(deg2rad(90 - XY(:, 2)), deg2rad(XY(:, 1)), ...
                    deg2rad(-lonc), deg2rad(latc), 0);
                lonp = rad2deg(phip);
                latp = 90 - rad2deg(thetap);
                [latf, lonf] = flatearthpoly(latp, lonp);
                XY = [lonf, latf];
            else

                if iscell(TH) && length(TH) >= 3
                    XY = feval(TH{:});
                else
                    XY = feval(dom, pars, buf, MoreRegionSpecs{:});
                end

            end

        case 'lonlat'
            % We have straight coordinates
            XY = TH;

        case 'polar'
            XY = [(1:360)', repmat(90 - TH, 360, 1)];
            [thetap, phip, ~] = ...
                rottp(deh2rad(90 - XY(:, 2)), deg2rad(XY(:, 1)), ...
                0, deg2rad(-theta), deg2rad(-phi));
            lonp = rad2deg(phip);
            latp = 90 - rad2deg(thetap);
            XY = [lonp latp];
    end

    % Run parallel is we are able
    parfor h = 1:J
        eigfunINT(h) = plm2avg(CC{h}, XY);
    end

    %% Returning requested outpput
    varargout = {eigfunINT};
end

%% Subfunctions
function varargout = parseinputs(inputArguments)
    CCDefault = '[~,CC]=localization(15,TH);';
    THDefault = 'africa';
    JDefault = [];
    phiDefault = 0;
    thetaDefault = 0;
    p = inputParser;
    addOptional(p, 'CC', CCDefault);
    addOptional(p, 'TH', THDefault, ...
        @(x) ischar(x) || isnumeric(x) || iscell(x) || isempty(x));
    addOptional(p, 'J', JDefault);
    addOptional(p, 'phi', phiDefault);
    addOptional(p, 'theta', thetaDefault);
    addOptional(p, 'MoreRegionSpecs', {});
    parse(p, inputArguments{:});

    CC = conddefval(p.Results.CC, CCDefault);
    TH = conddefval(p.Results.TH, THDefault);
    J = conddefval(p.Results.J, JDefault);
    phi = conddefval(p.Results.phi, phiDefault);
    theta = conddefval(p.Results.theta, thetaDefault);
    moreRegionSpecs = p.Results.MoreRegionSpecs;

    if iscell(TH)

        if length(TH) > 2
            moreRegionSpecs = {TH{3:end}, moreRegionSpecs{:}}; %#ok<CCAT>
        end

    end

    varargout = {CC, TH, J, phi, theta, moreRegionSpecs};
end

function regionType = regiontype(TH)

    if isnumeric(TH)

        if isscalar(TH)
            % We have a polar cap
            regionType = 'polar';
        else
            % We have straight coordinates
            regionType = 'lonlat';
        end

    elseif isstring(TH) || ischar(TH) || iscell(TH)
        % We have a geographic region
        regionType = 'geographic';
    else
        error('What format is TH?');
    end

end
