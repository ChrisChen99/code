function chrom = init(pop, jobNum, operaVec, candidateMachine, AGVNum, speedNum)
%% �� OS MS AS SS �Ĳ������
% OS���������
% MS����������
% AS��AGV����
% SS��AGV�ٶȵ�λ����
chrom = [];
operaNum = sum(operaVec);
opera = [];

%% OS�α���
for i = 1: jobNum
    opera = [opera, ones(1, operaVec(i)) * i];
end
% �� opera ���� pop ���Թ���һ�� pop*operaNum �ľ���
B = repmat(opera, pop, 1);
%����chaos
z_os = chaos(pop,operaNum);
% ��ÿһ�н�������
sorted_z = sort(z_os, 2);
% ��������һ�����󣬴洢����ĵȼ���
sorted_ranks = zeros(pop,operaNum);
for i = 1:pop
    [~, ranks] = ismember(z_os(i, :), sorted_z(i, :));  % �ҵ�ԭʼ������ÿ��Ԫ����������λ��
    sorted_ranks(i, :) = ranks;  % �洢����ĵȼ���
end
%%�洢��������OS����
OS = zeros(pop,operaNum);
for i = 1:pop
    [~, idx] = sort(sorted_ranks(i, :));
    OS(i, :) = B(i, idx);
end

%% MS�α���
z_ms = chaos(pop,operaNum);
%��ȡÿ������Ĺ����Ĺ���Ŀ��û����������Ͻ磩
MS = [];
upper_bound = [];
for i = 1: pop
    Mijk_num = [];
    for j = 1: jobNum
        for k = 1: operaVec(j)
            if isempty(candidateMachine{j, k})
                continue; % ������ڿ�ֵ��������ʣ�ಿ��
            end
            Mijk_num = [Mijk_num, length(candidateMachine{j, k})];
        end
    end
    upper_bound = repmat(Mijk_num, pop, 1);
end
MS = ceil(z_ms .* upper_bound);

%% AGV�����
z_as = chaos(pop,operaNum);
upper_bound = AGVNum;
AS = ceil(z_as .* upper_bound);

%% AGV�ٶȱ����
z_ss = chaos(pop,2*operaNum);
upper_bound = speedNum;
SS = ceil(z_ss .* upper_bound);

%% �ϲ�����5*operaNum��Ⱦɫ����Ⱥ
chrom = [OS,MS,AS,SS];
end

%% ����tentӳ�亯��
function z = chaos(pop, operaNum)
N = pop;
D = operaNum;
% 1) ���ɳ�ʼ��������� z1
z1 = rand(1, D);
% 2) ����ʣ�����
z = zeros(N, D);
z(1, :) = z1;
for i = 2:N
    for j = 1:D
        if z(i-1,j) < 0.5
            z(i,j) = 2 * z(i-1,j);
        else
            z(i,j) = 2 * (1 - z(i-1,j));
        end
        if z(i,j) == 0 || z(i,j) == 0.25 || z(i,j) == 0.5 || z(i,j) == 0.75
            z(i,j) = rand() * z(i,j);
        end
        if (i>1&&z(i,j)==z(i-1,j)) || (i>2&&z(i,j)==z(i-2,j)) ||...
                (i>3&&z(i,j)==z(i-3,j)) || (i>4&&z(i,j)==z(i-4,j))
            z(i,j) = rand() * z(i,j);
        end
    end
end
end