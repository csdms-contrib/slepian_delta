function varargout=filename2date(datanames,formatS,formatE)
% [thedates,datanames]=FILENAME2DATE(datanames,format)
%
% Extracts date information from filenames given format and location
% information.
%
%
% INPUT:
% 
% datanames    The name of the file you want to examine  OR
%              A cell array of filenames, as in ls2cell  OR
%              The full pathname to a directory with many files in it.  You
%               can use a regular expression (e.g. *) in the filename.
% formatS      A cell array where the first element is the character range 
%              in each filename that contains the date information 
%              (e.g. [17:22]).  The second element is the datenum format
%              string for this date (e.g. 'yyyymm').
% formatE      If you give a second format then it is assume that formatS
%              is a start date for the file and formatE is the end date for
%              the file.  This function will return thedates as the
%              midpoint of this timespan.
%
% OUTPUT:
%
% thedates     An array of the dates from the filenames in matlab datenum format
% datanames    The names of the netcdf files back to you.  Handy if you
%              gave a directory.
%    
%
% EXAMPLE: See default values
%
% Last modified by charig-at-princeton.edu, 09/02/2014

% Set some defaults as an example
defval('formatE',[]);
defval('datanames','GLDAS_NOAH10_M.A201401.001.nc');
defval('formatS',{[17:22] 'yyyymm'});


if iscell(datanames)
    % We have a list of files from ls2cell
    % This grabs the parts we want from all files at once
    %shortdata = cellfun(@(x) x(formatS{1}),datanames,'UniformOutput', false);
    %Then run them through datenum
    tempdateS = datenum(cellfun(@(x) x(formatS{1}),datanames,'UniformOutput', false)...
        ,formatS{2});
    
    if exist(formatE)
       tempdateE = datenum(cellfun(@(x) x(formatE{1}),datanames,...
           'UniformOutput', false),formatE{2});
    end
    
elseif exist(datanames)==2
    % It is one file, so just scan once
    tempdateS = datenum(datanames(formatS{1}),formatS{2});
    if exist(formatE)
        tempdateE = datenum(datanames(formatE{1}),formatE{2});
    end
    
elseif exist(datanames)==7
    % It is a directory, so we run ls2cell on the contents and then net2mat
    cls=ls2cell(datanames); 
    [tempdateS] = filename2date(cls,formatS,formatE);
    datanames = cls;
    
else
    error('Problem with your data name. Directories need fullfile paths!')
end

if exist(formatE)
    % Return the date as the midpoint between the beginning and enddates
    thedates = tempdateS' + (tempdateE' - tempdateS')/2;
else
    thedates = tempdateS';
end

    
% Collect output
varns={thedates,datanames};
varargout=varns(1:nargout);
    
   
    
    