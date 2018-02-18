function solution = WeightedAverage(data)

%This function will output the weighted average wave height of the given generalized
%solution U.

[weights,~] = randfixedsum(size(data.solution,2),1,1,0,1);

for i = 1:size(weights,1)
    weightedU(:,i) = weights(i).*data.solution(:,i);
end
solution = abs(sum(weightedU,2));