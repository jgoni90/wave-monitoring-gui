function [X,T,pos] = getSubmesh(Xd,Td)

%function [X,T,pos] = getSubmesh(Xd,Td)
%
%Returns the submesh X = Xd(pos,:), where pos = unique(Td), with
%connectivity T.

pos       = unique(Td);
X         = Xd(pos,:);
numv      = 1:max(pos);
numv(pos) = 1:size(X,1);
T         = numv(Td);