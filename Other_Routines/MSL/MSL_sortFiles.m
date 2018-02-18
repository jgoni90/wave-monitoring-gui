function [oFiles,n1,nend] = MSL_sortFiles(files)

%Argument "files" should be a cell of strings of the type "whatever052.mat"

nOfFiles = length(files);
pos = zeros(nOfFiles,1);
for i = 1:nOfFiles
    f = files{i};
    pos(i) = getFileNum(f);
end
[oPos,I] = sort(pos,'ascend');
n1 = oPos(1);
nend = oPos(end);
oFiles = files(I);

%%%%%%%%%%%%%%
function pos = getFileNum(file_name)

character = 1;
auxCounter = 4;
while ~isnan(character)
    character = str2double(file_name(end-auxCounter));
    auxCounter = auxCounter + 1;
end
pos = str2double(file_name(length(file_name)-auxCounter+2:end-4));
%%%%%%%%%%%%%%