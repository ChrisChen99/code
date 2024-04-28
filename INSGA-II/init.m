function chrom = init(pop, jobNum, operaVec, candidateMachine, AGVNum, speedNum)
%% 由 OS MS AS SS 四部分组成
% OS：工序编码
% MS：机器编码
% AS：AGV编码
% SS：AGV速度挡位编码
chrom = [];
operaNum = sum(operaVec);
opera = [];
for i = 1: jobNum
    opera = [opera, ones(1, operaVec(i)) * i];
end

for i = 1: pop
    OS = opera(randperm(operaNum));
    MS = [];
    for j = 1: jobNum
        for k = 1: operaVec(j)
            up = length(candidateMachine{j, k});
            MS = [MS, randperm(up, 1)];
        end
    end
    AS = randi([1, AGVNum], 1, operaNum);
    SS = randi([1, speedNum], 1, 2 * operaNum);
    chrom = [chrom; [OS, MS, AS, SS]];
end
end