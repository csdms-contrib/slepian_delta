function varargout = integratebasis(CC, TH, J, phi, theta)
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

    defval('TH', 'africa')
    defval('CC', '[~,CC]=localization(15,TH);');
    defval('phi', 0)
    defval('theta', 0)

    %%%
    % Initilization
    %%%

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
        if (floor(L) ~= L);
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
        [n, m] = size(CC);
        defval('J', m);
        defval('L', addmup(size(CC{1}, 1), 'r'));
    else
        error('What format is CC?');
    end

    % Now sort out what TH is
    if isstr(TH) % Geographic, we just do it
        % Lets check if we need to do a rotation. The function for your
        % coordinates should have this functionality if it's needed.
        defval('rotb', 0);

        try
            rotb = eval(sprintf('%s(''rotated'')', TH));
        end

        % Now do the rotation if needed
        if length(rotb) == 1 && rotb
            % Get the rotation parameters to rotate. Note, the region
            % rotation angles that we return from the functions (lonc, latc)
            % are the same regardless of if we did a buffer, as they pertain
            % to the original region
            [XY, lonc, latc] = eval(sprintf('%s(%i)', TH, 0));
            [thetap, phip, ~] = rottp((90 - XY(:, 2)) * pi / 180, XY(:, 1) * pi / 180, -lonc * pi / 180, latc * pi / 180, 0);
            lonp = phip * 180 / pi;
            latp = 90 - thetap * 180 / pi;
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
        defval('pars', 0);
        % Lets check if we need to do a rotation. The function for your
        % coordinates should have this functionality if it is needed.
        defval('rotb', 0);

        try
            rotb = eval(sprintf('%s(''rotated'')', dom));
        end

        % Now do the rotation if needed
        if length(rotb) == 1 && rotb
            % Return the coordinates and do the rotation
            [XY, lonc, latc] = eval(sprintf('%s(%i,%f)', TH{1}, 10, TH{2}));
            [thetap, phip, ~] = rottp((90 - XY(:, 2)) * pi / 180, XY(:, 1) * pi / 180, -lonc * pi / 180, latc * pi / 180, 0);
            lonp = phip * 180 / pi;
            latp = 90 - thetap * 180 / pi;
            [latf, lonf] = flatearthpoly(latp, lonp);
            XY = [lonf latf];
        else
            XY = eval(sprintf('%s(%i,%f)', dom, pars, buf));
        end

    elseif isfloat(TH) && isscalar(TH)
        % We have Polar caps
        XY = [(1:360)' repmat(90 - TH, 360, 1)];
        [thetap, phip, ~] = rottp((90 - XY(:, 2)) * pi / 180, XY(:, 1) * pi / 180, 0, -theta * pi / 180, -phi * pi / 180);
        lonp = phip * 180 / pi;
        latp = 90 - thetap * 180 / pi;
        XY = [lonp latp];
    else
        % Must be straight coordinates
        XY = TH;
    end

    % Initilization complete

    %%%
    % Do it
    %%%

    % Run parallel is we are able
    parfor h = 1:J
        [Int, ~] = plm2avg(CC{h}, XY);
        eigfunINT(h) = Int;
    end

    % Collect output
    varns = {eigfunINT};
    varargout = varns(1:nargout);
end
