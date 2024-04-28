function chrom = init(pop, jobNum, operaVec, candidateMachine, AGVNum, speedNum)
%% �� OS MS AS SS �Ĳ������
% OS���������
% MS����������
% AS��AGV����
% SS��AGV�ٶȵ�λ����
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