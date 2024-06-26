%% PopObj：算法求得的pareto解集
function Score = Spacing(PopObj)
Distance = pdist2(PopObj, PopObj, 'cityblock');
Distance(logical(eye(size(Distance, 1)))) = inf;
Score = std(min(Distance, [], 2));
end