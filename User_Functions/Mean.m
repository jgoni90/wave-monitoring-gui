function solution = Mean(data)

%This function will output the mean wave height of the given generalized
%solution U.

solution = abs(mean(data.solution,2));