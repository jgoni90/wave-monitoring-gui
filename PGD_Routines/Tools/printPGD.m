function printPGD(u,ue,i,ie)

if nargin == 0
    fprintf('\nMode\t\tResidual\tIterations\tPGD error\n')
    fprintf('----\t\t--------\t----------\t---------\n')
else
    fprintf('%i\t\t%e\t%i\t\t%e\n',u,ue,i,ie)
end