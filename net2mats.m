function varargout=net2mats(datanames,extension,varlist,savedir)
% [datanames,v]=NET2MATS(datanames,varlist,savedir)
%
% This program reads in netcdf files, extracts certain variables you want,
% returns the data, and saves them into a .mat file.  It can operate on a
% single file, multiple files, or a directory.
%
% This version is for the situation where you have many files
% and want to keep their .mat files separate. For example, to process 
% multiple GMT .grd files into .mat files. If you have many of the same
% files and want the same variable in each .mat file, then use NET2MAT (no S). 
% For example, to process a timeseries of GLDAS files to create a 
% timeseries of terrestrial water storage.
%
% INPUT:
% 
% datanames    The name of the netcdf file you want to examine.  Can also 
%              be a cell array of filenames, as in ls2cell.  Can also be 
%              the full pathname to a directory (ends in filesep) with 
%              many netcdf files in it.
% extension    The file extension as a string. [default: 'grd']
% varlist      A cell array of the variable list you want to extract and
%              save from the netcdf file(s). [default: All]
% savedir      Give alternate save location. Directories end with "filesep".
%              Default save location is in current working directory.
%
% OUTPUT:
%
% datanames    The names of the netcdf files back to you.  Handy if you
%              gave a directory.
% v            All of the variables you requested will be placed in a cell 
%              array in the second output. Each indicie in this cell array
%              is itself a cell array for each variable.  For example if
%              you read two netcdf files, requesting two variables, then
%              thevars will be 1x2, and thevars{1} will also be a 1x2 cell
%              array.
%
% EXAMPLE:
% [data1,data2] = net2mat('mynetcdffile.nc',{'variable1' 'variable2'});
%
% Last modified by charig-at-princeton.edu, 09/03/2014

defval('extension','grd');
defval('varlist',[]);
defval('savedir',[]);

% v will contain the variable names.  The netcdf variables themselves will be
% stored in variables of the same name, which are unkown to us.  Hence the
% use of eval statements on the contents of v.

if iscell(datanames)
    % A list of files
    % Each file goes in its own file, same as the name
    %v = genvarname(datanames);
    % Loop over the files
    for i=1:length(datanames)
        if isempty(varlist)
            % Use all the variables
            [varlist]=netvarread(datanames{i});
        end
        % If not empty, it should be a cell array
        %varnames = genvarname(varlist);

        % A list of variables, we need to loop over them
        for j=1:length(varlist)
            % Get data
            vardata = ncread(datanames{i},varlist{j});
            % Stick it in a variable with a good name
            eval([varlist{j} '= vardata;'])
        end

        % Now make a mat file
        [filepath,name,ext] = fileparts(datanames{i});
        save([filepath filesep name '.mat'],varlist{:})

    end
    disp(['NET2MAT: ' num2str(length(datanames)) ' data files read and'...
        ' saved as mat files.'])
        
elseif exist(datanames)==2
    % It is one file (given as a string), so just read it
    % ncread can only do one variable at a time, so we need to loop
    if isempty(varlist)
        % Use all the variables
        [varlist]=netvarread(datanames{i});
    end
    % If not empty, it should be a cell array
    %varnames = genvarname(varlist);

    % A list of variables, we need to loop over them
    for j=1:length(varlist)
        % Get data
        vardata = ncread(datanames{i},varlist{j});
        % Stick it in a variable with a good name
        eval([varlist{j} '= vardata;'])
    end

    % Now make a mat file
    [filepath,name,ext] = fileparts(datanames{i});
    save([filepath filesep name '.mat'],varlist{:})

elseif exist(datanames)==7
    % It is a directory, so we run ls2cell on the contents and then net2mats
    defval('savedir',datanames);
    cls=ls2cell([datanames '*' extension],1); 
    % We want only certain files in there.
    [datanames,thevars] = net2mats(cls,varlist,savedir);

    %for k=1:length(varlist); eval([varlist{k} '= thevars{k};']); end;
    
    
else
    error('Problem with your data name. Directories need fullfile paths!')
end

% % Create the mat files
varns{1} = datanames;
% for j=1:length(varlist)
%     if exist(savedir)
%         % Directories already end with "filesep"
%         save([savedir v{j} '.mat'],v{j}); 
%     else
%         save([v{j} '.mat'],v{j});
%     end
%     varns{j+1} = eval([v{j}]); 
% end

    
% Provide all the output
varargout{1}=varns{1};
varargout{2} = {varns{2:end}};



   
    
    