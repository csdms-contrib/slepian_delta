function varargout = gracedeg2(Rlevel)
    % [deg2data]=GRACEDEG2(Rlevel)
    %
    % This function reads and formats the degree-2 spherical harmonic
    % correction from GRACE/GRACE-FO Technical Note 14.
    %
    % INPUT:
    %
    % Rlevel      The release level of the solution you want.
    %              Either 'RL04','RL05', or 'RL06'
    %              Note: currently only RL06 supported
    %
    % OUTPUT:
    %
    % Returns these variables and saves them in a .mat file:
    % deg2data       Three column matrix where the first column is dates, the
    %                second column is C20 values, and the third column is
    %                C30 values.
    %
    % NOTE:
    %
    % The GRACE degree-2 correction data are distributed in Technical Note 14,
    %   TN-14. You can find TN-14 in the docs folder in the GRACE data directory at NASA
    %   PODAAC. https://podaac.jpl.nasa.gov/
    %
    % We check the file timestamps of both the TN-14 file and the Matlab .mat
    %   file that we save/load. If the .mat file is newer, then we load it. If
    %   the TN-14 document is newer, then we use it to create a new .mat file
    %   and then load the .mat file.
    %
    % These
    %
    % Last modified by charig-at-arizona.edu, 11/11/2021

    % Determine parameters and set defaults
    defval('Rlevel', 'RL06')

    % Where the original data files are kept
    defval('ddir1', fullfile(getenv('IFILES'), 'GRACE', 'Degree2'));
    % Where you would like to save the new .mat file
    defval('ddir2', fullfile(getenv('IFILES'), 'GRACE', 'Degree2'));

    % The name of the original data files
    if strcmp(Rlevel, 'RL06')
        fnpl1 = sprintf('%s/TN-14_C30_C20_GSFC_SLR.txt', ddir1);
        fnpl2 = sprintf('%s/TN-14_C30_C20_GSFC_SLR_%s.mat', ddir2, Rlevel);
    elseif strcmp(Rlevel, 'RL05')
        %fnpl1=sprintf('%s/deg1_RL05_NH.txt',ddir1);
        %fnpl2=sprintf('%s/deg1_RL05_NH.mat',ddir2);
    elseif strcmp(Rlevel, 'RL04')
        %fnpl1=sprintf('%s/deg1_RL04_NH.txt',ddir1);
        %fnpl2=sprintf('%s/deg1_RL04_NH.mat',ddir2);
    else
        error('GRACEDEG2: Wonky release level requested')
    end

    % Get the file date information
    d1 = dir(fnpl1);
    d2 = dir(fnpl2);

    % If this file already exists, and it is more new than the source file,
    % then load it.  Otherwise make a new one and return that one. This should
    % automatically make and use a new .mat file if you update the source (TN-14) file.
    if ~(exist(fnpl2, 'file') == 2) || datenum(d1.date) > datenum(d2.date)
        % Here create a new save file to use.
        % This is reversed from our typical if, exist, load convention because
        % this way if the file does not exist it will skip the second if
        % condition (which would error on the datenum call) and go straight here

        fid = fopen(fnpl1);
        tline = fgetl(fid);

        % Keep reading lines until we get to the 'Product:' line
        while isempty(tline) || ~strcmp(tline(1:6), 'Produc')
            % Just 6 letters here because line 16 of this file is "Notes:" and
            % comparing 7 would error.
            tline = fgetl(fid);
        end

        % Now do a formatted text scan on the data (10 columns)
        C = textscan(fid, '%f%f%f%f%f%f%f%f%f%f');
        fclose(fid);

        % Use the epoch start and end dates to get the midpoint date
        startepoch = datenum(datetime(C{1}, 'ConvertFrom', 'modifiedjuliandate'));
        endepoch = datenum(datetime(C{9}, 'ConvertFrom', 'modifiedjuliandate'));

        % Make sure we have unique dates
        thedates = (startepoch + endepoch) ./ 2;
        [b, m] = unique(thedates);
        thedates = thedates(m);

        %     % Now reorder into a lmcosi format for each month
        %      for i=1:length(thedates)
        %          deg2lmcosi(i,:,:) = [2 0 C{3}(i) 0; 3 0 C{6}(i) 0];
        %      end

        deg2data = [thedates, C{3}, C{6}];

        % Create a save file
        save(fnpl2, 'deg2data')

    else
        load(fnpl2)
        fprintf('%s loading %s\n', upper(mfilename), fnpl2)
    end

    % Collect output
    varns = {deg2data};
    varargout = varns(1:nargout);
end