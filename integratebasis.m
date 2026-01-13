%% INTEGRATEBASIS
% [eigfunINT]=INTEGRATEBASIS(CC,TH,J)
%
% Accepts a Slepian basis and integrates the functions within a region.
% The integrals reported are in fractional sphere area.  You may want to
% convert this to real sphere area via x(4*pi*radius^2).
% For geographic regions, the function will automatically save the data
% to a file in the IFILES directory.
%
% INPUT:
%
% CC        The eigenfunctions.  This can be either the G you get from
%           GLMALPHA or a cell array of the functions as in LOCALIZATION.
% TH        The region.  Can be formatted multiple ways as in GLMALPHA
% J         How many functions you want to do this for [DEFAULT: all]
% phi       Longitude of the center (degrees)
% theta     Colatitude of the center (degrees)
%           Note: there is no omega here because a rotation of the
%           polar cap functions does not affect their integration
% ForceNew  Force the function to recalculate the integrals [DEFAULT: false]
% SaveData  Save the data to a file [DEFAULT: true]
%
% OUTPUT:
%
% eigfunINT The integrals of the eigenfunctions over the region.  These
%           are in real sphere area, depending on your radius.
%
% SEE ALSO: PLM2AVG, PLM2AVGP
%
% Last modified by
%   williameclee-at-arizona.edu, 10/31/2024
%   charig-at-email.arizona.edu, 11/01/2016

function eigfunINT = integratebasis(varargin)
    [CC, TH, J, phi, theta, forceNew, saveData] = parseinputs(varargin{:});
    defval('TH', 'africa')
    defval('CC', '[~,CC]=localization(15,TH);');
    defval('phi', 0)
    defval('theta', 0)

    %% Initilisation
    % Sort out what CC is
    if isstring(CC) || ischar(CC)
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
        if floor(L) ~= L
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

    % Now sort out what TH is
    if isstring(TH) || ischar(TH) % Geographic, we just do it
        % Lets check if we need to do a rotation. The function for your
        % coordinates should have this functionality if it's needed.
        defval('rotb', 0);

        try
            rotb = feval(TH, 'rotated');
        catch
        end

        % Now do the rotation if needed
        if isscalar(rotb) && rotb
            % Get the rotation parameters to rotate. Note, the region
            % rotation angles that we return from the functions (lonc, latc)
            % are the same regardless of if we did a buffer, as they pertain
            % to the original region
            [XY, lonc, latc] = feval(TH, 10);
            [thetap, phip, ~] = rottp(deg2rad(90 - XY(:, 2)), deg2rad(XY(:, 1)), deg2rad(-lonc), deg2rad(latc), 0);
            lonp = rad2deg(phip);
            latp = 90 - rad2deg(thetap);
            [latf, lonf] = flatearthpoly(latp, lonp);
            XY = [lonf latf];
        else
            % No changes
            XY = TH;
        end

    elseif iscell(TH) % Geographic + buffer
        defval('buf', 0);
        dom = TH{1};
        buf = TH{2};
        defval('pars', 10);
        % Lets check if we need to do a rotation. The function for your
        % coordinates should have this functionality if it is needed.
        defval('rotb', 0);

        try
            rotb = feval(dom, 'rotated');
        catch
        end

        % Now do the rotation if needed
        if isscalar(rotb) && rotb
            % Return the coordinates and do the rotation
            [XY, lonc, latc] = feval(TH{1}, 10, TH{2});
            [thetap, phip, ~] = rottp(deg2rad(90 - XY(:, 2)), deg2rad(XY(:, 1)), deg2rad(-lonc), deg2rad(latc), 0);
            lonp = rad2deg(phip);
            latp = 90 - rad2deg(thetap);
            [latf, lonf] = flatearthpoly(latp, lonp);
            XY = [lonf latf];
        else
            XY = feval(dom, pars, buf);
        end

    elseif isfloat(TH) && isscalar(TH)
        % We have Polar caps
        XY = [(1:360)' repmat(90 - TH, 360, 1)];
        [thetap, phip, ~] = rottp(deg2rad(90 - XY(:, 2)), deg2rad(XY(:, 1)), deg2rad(-theta), deg2rad(phi), 0);
        lonp = rad2deg(phip);
        latp = 90 - rad2deg(thetap);
        XY = [lonp latp];
    else
        % Must be straight coordinates
        XY = TH;
    end

    % Initilization complete

    %% Find if we have saved data
    [dataPath, hasDataSaved] = datapath(TH, J, phi, theta);

    if hasDataSaved && ~forceNew
        load(dataPath, 'eigfunINT');
        fprintf('%s loaded %s\n', upper(mfilename), dataPath);
        return
    end

    %% Do it
    % Run parallel is we are able
    parfor h = 1:J
        eigfunINT(h) = plm2avg(CC{h}, XY);
    end

    %% Save
    if dataPath && saveData
        save(dataPath, 'eigfunINT');
        fprintf('%s saved %s\n', upper(mfilename), dataPath);
    end

end

%% Subfunctions
function varargout = parseinputs(varargin)
    ip = inputParser;
    addOptional(ip, 'CC', [], @(x) isnumeric(x) || iscell(x) || ischar(x) || isempty(x));
    addOptional(ip, 'TH', [], @(x) isnumeric(x) || iscell(x) || ischar(x) || isempty(x));
    addOptional(ip, 'J', [], @(x) isnumeric(x) || isempty(x));
    addOptional(ip, 'phi', [], @(x) isnumeric(x) || isempty(x));
    addOptional(ip, 'theta', [], @(x) isnumeric(x) || isempty(x));
    addParameter(ip, "ForceNew", false, @(x) islogical(x) || isnumeric(x));
    addParameter(ip, "SaveData", true, @(x) islogical(x) || isnumeric(x));
    parse(ip, varargin{:});
    CC = ip.Results.CC;
    TH = ip.Results.TH;
    J = ip.Results.J;
    phi = ip.Results.phi;
    theta = ip.Results.theta;
    forceNew = ip.Results.ForceNew;
    saveData = ip.Results.SaveData;
    varargout = {CC, TH, J, phi, theta, forceNew, saveData};
end

function [dataPath, hasDataSaved] = datapath(TH, J, phi, theta)

    if ~ismember(class(TH), {'cell', 'string', 'char'}) || (isnumeric(TH) && ~isscalar(TH))
        dataPath = false;
        hasDataSaved = false;
        return
    end

    dataFolder = fullfile(getenv('IFILES'), 'EIGFUNINT');

    if ~exist(dataFolder, 'dir')
        mkdir(dataFolder);
        fprintf('%s created directory %s\n', upper(mfilename), dataFolder);
    end

    if iscell(TH) && TH{2} ~= 0
        dataFile = sprintf('EIGFUNINT-%s-%f-%d.mat', TH{1}, TH{2}, J);
    elseif iscell(TH)
        TH = TH{1};
        dataFile = sprintf('EIGFUNINT-%s-%d.mat', TH, J);
    elseif ischar(TH) || isstring(TH)
        dataFile = sprintf('EIGFUNINT-%s-%d.mat', TH, J);
    elseif isnumeric(TH)
        dataFile = sprintf('EIGFUNINT-%f-%f%f-%d.mat', TH, phi, theta, J);
    end

    dataPath = fullfile(dataFolder, dataFile);
    hasDataSaved = exist(dataPath, 'file') == 2;
end
