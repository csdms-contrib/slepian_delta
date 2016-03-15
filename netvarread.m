function varargout=netvarread(myfile)
% [varlist,allsame]=NETVARREAD(myfile)
%
% Read a netcef file and get the variable list in the file.
%
%
% INPUT:
% 
% myfile       The name of file you want to examine. OR  
%              A list of filenames.  OR
%              A directory containing netcdf files.
%                If all have the same variables you get one list.
%
% OUTPUT:
%     
% varlist      The variable name list in a cell array
% allsame      If you gave a cell array of files, did they all have the
%                same variables?
%
%
% Last modified by charig-at-princeton.edu, 09/02/2014

defval('varlist',[]);
defval('allsame',1);


if iscell(myfile)
    % We have a list of files, e.g. from ls2cell
    Finfo = ncinfo(myfile{1});
    filevarlist{1} = {Finfo.Variables.Name};
    
    for i = 2:length(myfile)
        Finfo = ncinfo(myfile{i});
        filevarlist{i} = {Finfo.Variables.Name};
        if ~isequal(filevarlist{1},filevarlist{i})
            allsame=0;
        end
    end
    
    if allsame
        %disp(['All of your files have the same variable names in them. Only'... 
        %' one list returned.'])
        varlist = filevarlist{1};
    else
        varlist = filevarlist;
    end
   
elseif exist(myfile)==2
    % It is one file, so just read it
    %ischar(myfile)
    Finfo = ncinfo(myfile);
    varlist = {Finfo.Variables.Name};
    
elseif exist(myfile)==7
    % It is a directory, so we run ls2cell on the contents and then netvarread
    
    cls=ls2cell([myfile '*nc'],1); % We want only the netcdf files in there.
    [varlist,allsame] = netvarread(cls);
    
end

    
% Provide output where requested
varns={varlist,allsame};
varargout=varns(1:nargout);
    
   
    
    