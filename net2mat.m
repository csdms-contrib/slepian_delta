function varargout=net2mat(datanames,varlist,savedir)
% [datanames,v]=NET2MAT(datanames,varlist,savedir)
%
% This program reads in netcdf files, extracts certain variables you want,
% returns the data, and saves them into a .mat file.  It can operate on a
% single file, multiple files, or a directory.
%
%
% INPUT:
% 
% datanames    The name of the netcdf file you want to examine.  Can also 
%              be a cell array of filenames, as in ls2cell.  Can also be 
%              the full pathname to a directory (ends in filesep) with 
%              many netcdf files in it.
% varlist      A cell array of the variable list you want to extract and
%              save from the netcdf file(s).
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

defval('varlist',[]);
defval('allsame',1);
defval('savedir',[]);

% v will contain the variable names.  The netcdf variables themselves will be
% stored in variables of the same name, which are unkown to us.  Hence the
% use of eval statements on the contents of v.

if iscell(datanames)
    % Each variable goes in its own file, same as the name
    v = genvarname(varlist);
    % We have a list of files from ls2cell, need more loops
    for i=1:length(datanames)
        if ischar(varlist)
           % One variable in each file will be put into one variable, named
           % same as the varlist, which is a cell array with an indicie for
           % each file.
           vardata = ncread(datanames{i},varlist);
           eval([v '{i}= vardata;']);
        elseif iscell(varlist)
           % A list of variables, we need to loop over them
           for j=1:length(varlist)
               % Same as above but make multiple variables from varlist
               vardata = ncread(datanames{i},varlist{j});
               eval([v{j} '{i}= vardata;']);
           end
        end
    end
    disp(['NET2MAT: ' num2str(length(datanames)) ' data files read and'...
        ' saved as variable files.'])
        
elseif exist(datanames)==2
    % It is one file, so just read it
    % ncread can only do one variable at a time, so we need to loop
    if ischar(varlist)
        % One variable, just read it
        vardata = ncread(datanames,varlist);
    elseif iscell(varlist)
        % A list of variables, we need to loop over them
        % Each variable goes in its own file, same as the name
        v = genvarname(varlist);
        for j=1:length(varlist)
            vardata = ncread(datanames,varlist{j});
            eval([v{j} '= vardata;']);
        end       
    end

elseif exist(datanames)==7
    % It is a directory, so we run ls2cell on the contents and then net2mat
    defval('savedir',datanames);
    v = genvarname(varlist);
    %command1 = ['[' varlist{1}];
    %if length(varlist)>1
    %   for k=2:length(varlist); command1 = [command1 ',' varlist{k}]; end;
    %end
    %command1 = [command1 ']'];
    cls=ls2cell([datanames '*nc'],1); % We want only the netcdf files in there.
    %eval([command1 '=net2mat(cls,varlist)']);
    [datanames,thevars] = net2mat(cls,varlist,savedir);
    %datanames = cls;
    for k=1:length(varlist); eval([varlist{k} '= thevars{k};']); end;
    
    
else
    error('Problem with your data name. Directories need fullfile paths!')
end

% Create the mat files
varns{1} = datanames;
for j=1:length(varlist)
    if exist(savedir)
        % Directories already end with "filesep"
        save([savedir v{j} '.mat'],v{j}); 
    else
        save([v{j} '.mat'],v{j});
    end
    varns{j+1} = eval([v{j}]); 
end

    
% Provide all the output
varargout{1}=varns{1};
varargout{2} = {varns{2:end}};



   
    
    