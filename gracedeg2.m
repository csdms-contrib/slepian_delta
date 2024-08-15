%% GRACEDEG2
% Reads and formats the C20 and C30 spherical harmonic corrections from GRACE/GRACE-FO Technical Note 14.
%
% Syntac
%   c2030data = gracedeg2(Rlevel)
%   [c2030data, c2030sigma] = gracedeg2(Rlevel)
%
% Inputs
%   Rlevel - The release level of the solution you want.
%       Either 'RL04','RL05', or 'RL06'
%       Note: currently only RL06 supported
%
% Outputs and saved variables
%   c2030data - The C20 and C30 coefficients from SLR data
%       The first column is the time (in datenum)
%       The second column is the C20 coefficient
%       The third column is the C30 coefficient (only for GRACE-FO)
%   c2030sigma - The C20 and C30 sigma from SLR data
%       The first column is the time (in datenum)
%       The second column is the C20 sigma
%       The third column is the C30 sigma (only for GRACE-FO)
%
% Note
%   The GRACE degree-2 correction data are distributed in Technical Note 14,
%   TN-14. You can find TN-14 in the docs folder in the GRACE data directory at NASA
%   PODAAC. https://podaac.jpl.nasa.gov/
%   We check the file timestamps of both the TN-14 file and the Matlab .mat
%   file that we save/load. If the .mat file is newer, then we load it. If
%   the TN-14 document is newer, then we use it to create a new .mat file
%   and then load the .mat file.
%
% Last modified by 
%   williameclee-at-arizona.edu, 07/30/2024
%   charig-at-arizona.edu, 11/11/2021

function varargout = gracedeg2(varargin)
    p = inputParser;
    p.KeepUnmatched = true;
    addOptional(p, 'Rlevel', 'RL06', ...
        @(x) ischar(x) && ismember(x, {'RL04', 'RL05', 'RL06'}));
    addParameter(p, 'ForceNew', false, @islogical);
    parse(p, varargin{:});
    Rlevel = p.Results.Rlevel;
    forceNew = p.Results.ForceNew;

    % Where the original data files are kept
    inputFolder = fullfile(getenv('IFILES'), 'GRACE', 'Degree2');
    % Where you would like to save the new .mat file
    outputFolder = fullfile(getenv('IFILES'), 'GRACE', 'Degree2');

    % The name of the original data files
    switch Rlevel
        case 'RL06'
            inputPath = fullfile(inputFolder, ...
            'TN-14_C30_C20_GSFC_SLR.txt');
            outputPath = fullfile(outputFolder, ...
                sprintf('TN-14_C30_C20_GSFC_SLR_%s.mat', Rlevel));
        case 'RL05'
            inputPath = fullfile(inputFolder, 'deg1_RL05_NH.txt');
            outputPath = fullfile(outputFolder, 'deg1_RL05_NH.mat');
        case 'RL04'
            inputPath = fullfile(inputFolder, 'deg1_RL04_NH.txt');
            outputPath = fullfile(outputFolder, 'deg1_RL04_NH.mat');
        otherwise
            error('Wonky release level requested')
    end

    % Get the file date information
    d1 = dir(inputPath);
    d2 = dir(outputPath);

    % If this file already exists, and it is more new than the source file,
    % then load it. Otherwise make a new one and return that one. This 
    % should automatically make and use a new .mat file if you update the 
    % source (TN-14) file.
    if exist(outputPath, 'file') == 2 && d2.datenum > d1.datenum ...
            && ~forceNew
        fprintf('%s loading %s\n', upper(mfilename), outputPath)
        load(outputPath, 'c2030data', 'c2030sigma')

        if exist('c2030data', 'var') && exist('c2030sigma', 'var')
            varargout = {c2030data, c2030sigma};
            return
        end

    end

    % Here create a new save file to use.
    % This is reversed from our typical if, exist, load convention because
    % this way if the file does not exist it will skip the second if
    % condition (which would error on the datenum call) and go straight here

    fid = fopen(inputPath);
    tline = fgetl(fid);

    % Keep reading lines until we get to the 'Product:' line
    while ~contains(tline, 'Product:')
        % Just 6 letters here because line 16 of this file is "Notes:" and
        % comparing 7 would error.
        tline = fgetl(fid);
    end

    % Now do a formatted text scan on the data (10 columns)
    C = textscan(fid, '%f%f%f%f%f%f%f%f%f%f');
    fclose(fid);

    % Use the epoch start and end dates to get the midpoint date
    startdates = datenum(datetime(C{1}, ...
        'ConvertFrom', 'modifiedjuliandate')); %#ok<DATNM>
    enddates = datenum(datetime(C{9}, ...
        'ConvertFrom', 'modifiedjuliandate')); %#ok<DATNM>

    % Make sure we have unique dates
    dates = (startdates + enddates) ./ 2;
    [~, m] = unique(dates);
    dates = dates(m);
    c2030data = [dates, C{3}, C{6}];
    c2030sigma = [dates, C{5} * 1e-10, C{8} * 1e-10];

    % Create a save file
    save(outputPath, 'c2030data', 'c2030sigma')

    % Collect output
    varargout = {c2030data, c2030sigma};
end
