% [thedates,deg1data]=GRACEDEG1(Pcenter,Rlevel)
%
% This function reads and formats the degree-1 spherical harmonic
% correction from GRACE/GRACE-FO Technical Note 13.
%
% INPUT:
%
% Pcenter     'CSR' data center at the Center for Space Research
%             'GFZ' data center at the GeoForschungsZentrum Potsdam
%             'JPL' data center at the Jet Propulsion Laboratory
% Rlevel      The release level of the solution you want.
%              Either 'RL04','RL05', or 'RL06'
%
% OUTPUT:
%
% Returns these variables and saves them in a .mat file:
%    thedates       time stamps in Matlab time. This date is the midpoint
%                      of the GRACE month data span.
%    deg1data       potential coefficients for just degree l = 1
%                      [nmonths x 2 x 4]
%
% NOTES:
%
% The GRACE degree-1 correction data are distributed in Technical Note 13,
%   TN-13. There is one for each datacenter because the GRACE coefficients
%   from degree 2-60 are used in the construction of degree one terms. You
%   can find TN-13 in the docs folder in the GRACE data directory at NASA
%   PODAAC. https://podaac.jpl.nasa.gov/
%
% We check the file timestamps of both the TN-13 file and the Matlab .mat
%   file that we save/load. If the .mat file is newer, then we load it. If
%   the TN-13 document is newer, then we use it to create a new .mat file
%   and then load the .mat file.
%
% Last modified by charig-at-arizona.edu, 11/18/2021

function varargout = gracedeg1(varargin)
    p = inputParser;
    p.KeepUnmatched = true;
    addOptional(p, 'Pcenter', 'CSR', ...
        @(x) ischar(x) && ismember(x, {'CSR', 'GFZ', 'JPL'}));
    addOptional(p, 'Rlevel', 'RL06', ...
        @(x) ischar(x) && ismember(x, {'RL04', 'RL05', 'RL06'}));
    addParameter(p, 'ForceNew', false, @islogical);

    parse(p, varargin{:});
    Pcenter = p.Results.Pcenter;
    Rlevel = p.Results.Rlevel;
    forceNew = p.Results.ForceNew;

    % Determine parameters and set defaults
    defval('Pcenter', 'CSR')
    defval('Rlevel', 'RL06')

    % Where the original data files are kept
    inputFolder = fullfile(getenv('IFILES'), 'GRACE', 'Degree1');
    % Where you would like to save the new .mat file
    outputFolder = fullfile(getenv('IFILES'), 'GRACE', 'Degree1');

    % The name of the original data files
    switch Rlevel
        case 'RL06'
            inputPath = fullfile(inputFolder, ...
                sprintf('TN-13_GEOC_%s_%s.txt', Pcenter, 'RL0602'));
            outputPath = fullfile(outputFolder, ...
                sprintf('TN-13_GEOC_%s_%s.mat', Pcenter, 'RL0602'));
        case 'RL05'
            inputPath = fullfile(inputFolder, ...
                sprintf('deg1_%s_NH.txt', 'RL05'));
            outputPath = fullfile(outputFolder, ...
                sprintf('deg1_%s_NH.mat', 'RL05'));
        case 'RL04'
            inputPath = fullfile(inputFolder, ...
                sprintf('deg1_%s_NH.txt', 'RL04'));
            outputPath = fullfile(outputFolder, ...
                sprintf('deg1_%s_NH.mat', 'RL04'));
        otherwise
            error('GRACEDEG1: Wonky release level requested')
    end

    % Get the file date information
    d1 = dir(inputPath);
    d2 = dir(outputPath);

    % If this file already exists, and it is more new than the source file,
    % then load it.  Otherwise make a new one and return that one. This should
    % automatically make and use a new .mat file if you update the source (TN-13) file.
    if exist(outputPath, 'file') == 2 && d2.datenum > ...
            d1.datenum && ~forceNew
        fprintf('%s loading %s\n', upper(mfilename), outputPath)
        load(outputPath, 'thedates', 'deg1data')

        varargout = {thedates, deg1data};
        return
    end

    % Here create a new save file to use.
    % This is reversed from our typical if, exist, load convention because
    % this way if the file does not exist it will skip the second if
    % condition (which would error on the datenum call) and go straight here

    fid = fopen(inputPath);
    tline = fgetl(fid);

    % Keep reading lines until we get to the 'end of header' line
    finalHeaderLine = 'end of header';

    while ~contains(tline, finalHeaderLine)
        tline = fgetl(fid);
    end

    clear d1 d2

    %% Reading the data
    % Now do a formatted text scan on the data
    C = textscan(fid, '%s%d%d%f%f%f%f%s%s%d');
    % Use string format for the dates so we can easily send it to datenum
    fclose(fid);

    % What we actually want is a 3 dimensional cell array with the date as
    % the first dimension and a lmcosi as dimension 2 and 3 (just like
    % grace2plmt)

    % Use the epoch start and end dates to get the midpoint date
    thedates = (datenum(C{8}, 'yyyymmdd') + datenum(C{9}, 'yyyymmdd')) ./ 2; %#ok<DATNM>
    [~, m] = unique(thedates);
    thedates = thedates(m);

    % Re cast the ints as doubles to get it all in one matrix
    deg1data = [double(C{2}), double(C{3}), C{4}, C{5}];
    deg1lmcosi = zeros(length(C{1}) / 2, 2, 4);

    % Now reorder into a lmcosi format for each month
    for i = 1:(length(C{1})) / 2
        temp = [deg1data(2 * i - 1, :);
                deg1data(2 * i, :)];
        deg1lmcosi(i, :, :) = temp;
    end

    deg1data = deg1lmcosi;

    % Create a save file
    save(outputPath, 'thedates', 'deg1data')

    % Collect output
    varargout = {thedates, deg1data};
end
