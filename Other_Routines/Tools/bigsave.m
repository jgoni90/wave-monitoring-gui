function bigsave(file,varargin)

% bigsave(file,'var1','var2',...) saves the variables var1,
% var2,... independently into the mat files file_var1.mat,file_var2.mat,...
% using the matlab flag -v7.3 for those varibles with size >= 2Gb.
%
% bigsave(file,'-struct','s') saves the fields of structure s as
% independent variables into mat files file_field1.mat,file_field2.mat,...
% using the matlab flag -v7.3 for those varibles with size >= 2Gb.
%
% June 2013 by David Modesto (david.modesto@upc.edu, LaCàN, UPC)

%Check input
if (nargin == 3) && strcmpi(varargin{1},'-struct')
    structflag = true;
else
    structflag = false;
end

%Max gb allowed by Matlab without using -v7.3 flag
maxgb = 1.99;

%Save
if structflag %for variables which are structure fields
    structvar = evalin('caller',varargin{2});
    structnames = fieldnames(structvar);
    nOfsaves = numel(structnames);
    for i = 1:nOfsaves
        ifile = [file '_' structnames{i} '.mat'];
        eval([structnames{i} '= structvar.(structnames{i});']);
        gb = getGb(structvar.(structnames{i}));
        if gb >= maxgb
            disp(['Saving variable ' structnames{i} ' with size ' num2str(gb)...
                ' using flag -v7.3...'])
            save(ifile,structnames{i},'-v7.3')
        else
            save(ifile,structnames{i})
        end
    end
else %for independent variables
    nOfsaves = numel(varargin);
    for i = 1:nOfsaves
        if isempty(varargin{i}), continue, end
        ifile = [file '_' varargin{i} '.mat'];
        ivar = evalin('caller',varargin{i});
        eval([varargin{i} '= ivar;']);
        gb = getGb(ivar);
        if gb >= maxgb
            disp(['Saving variable ' varargin{i} ' with size ' num2str(gb)...
                ' using flag -v7.3...'])
            save(ifile,varargin{i},'-v7.3')
        else
            save(ifile,varargin{i})
        end
    end
end

%Function to get size in Gb 
function gb = getGb(var)

varinfo = whos('var');
gb = varinfo.bytes / (1024^3);
