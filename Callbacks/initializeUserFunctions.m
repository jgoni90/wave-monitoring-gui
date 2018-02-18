function initializeUserFunctions(handles)

directory = [pwd '/User_Functions'];
functions = dir([directory '/*.m']);
[sorted_names,~] = sortrows({functions.name}');
functionnames = regexprep(sorted_names,'.m','');
set(handles.UserFcn,'String',functionnames)