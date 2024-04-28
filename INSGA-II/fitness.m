function [FUNC, machineTable, AGVTable, makespan, EG_M_SUM, EG_A_SUM, agvEGRecord, agvChargeNum] = ...
    fitness(chrom, jobNum, jobInfo, operaVec, machineNum, AGVNum, AGVSpeed, candidateMachine, ...
    distance_matrix, machineEnergy, AGVEnergy, AGVEG_MAX, AGVEG_MIN, eChargeSpeed)
%% machineTable ��¼ÿ�������Ĺ���ʱ���
% �����ṹ�� workTable
workTable = [];
workTable.start = 0;    % ��ʼʱ��
workTable.end = Inf;    % ����ʱ��
workTable.job = 0;      % �ӹ��Ĺ����ţ�0��ʾ����
workTable.opera = 0;    % �ӹ��Ĺ����
% ��ʼ��machineTable
for i = 1 : machineNum
    machineTable{1, i} = workTable;
end

%% AGVTable ��¼ÿ��AGV���ƶ�ʱ���
% �����ṹ�� transferTable
transferTable = [];
transferTable.start = 0;
transferTable.end = Inf;
transferTable.job = 0;          % ���˵Ĺ�����
transferTable.opera = 0;        % ���˵Ĺ����
transferTable.load_status = 0;  % ���˵�״̬   -1 ����ת��  -2 ����ת��
transferTable.from_machine = -1;    % ���ĸ����������� -1��װ��վ  -2��ж��վ  -3�����׮��ж��վ�� 0������ �������������
transferTable.to_machine = 0;       % ͬ��
transferTable.charge = 0;       % 0 ����״̬  1 ���״̬  2 ǰ�����״̬
% ��ʼ�� AGVTable
for i = 1 : AGVNum
    AGVTable{1, i} = transferTable;
end

%% ��������
[machineTable, AGVTable, jobCompleteUnLoad, agvEGRecord, agvChargeNum] = sorting(chrom, jobNum, jobInfo, operaVec, AGVNum, ...
    AGVSpeed, candidateMachine, distance_matrix, AGVEnergy, AGVEG_MAX, AGVEG_MIN, eChargeSpeed, machineTable, AGVTable);

% �������Լ��ʹ��
% figure(1)
% machine_AGV_gantt_chart(machineTable, AGVTable);
% figure(2)
% ag_str = [];
% for ag = 1: AGVNum
%     ag_str = [ag_str ',''AGV' num2str(ag), ''''];
% end
% ag_str = ag_str(2: end);
%
% color = ['r'; 'g'; 'b'; 'y'; 'c'];
% for ag = 1: AGVNum
%     plot(agvEGRecord{ag}(:, 1), agvEGRecord{ag}(:, 2), [color(ag) '-s'], 'MarkerEdgeColor', 'k', 'MarkerFaceColor', color(ag), 'LineWidth', 1.0);
%     hold on
%     for i = 1: size(agvEGRecord{ag}, 1)
%         txt = sprintf('%0.2f Kw', agvEGRecord{ag}(i, 2));
%         text(agvEGRecord{ag}(i, 1) + 1, agvEGRecord{ag}(i, 2) - 1, txt, 'FontWeight', 'Bold', 'FontSize', 8)
%         hold on
%     end
% end
% eval(['legend(' ag_str ', ''Location'', ''NorthEastOutside'')'])

%% Ŀ�꺯�� 1 ��makespan
makespan = max(jobCompleteUnLoad);

%% Ŀ�꺯�� 2�����ܺ�
%
% [1]. �����ܺ�
machine_work = zeros(machineNum, 1);
machine_spare = zeros(machineNum, 1);
for i = 1: machineNum
    for j = 1: length(machineTable{i})
        if isequal(machineTable{i}(j).end, inf)
            continue;
        end

        if isequal(machineTable{i}(j).job, 0)
            machine_spare(i) = machine_spare(i) + (machineTable{i}(j).end - machineTable{i}(j).start);
        else
            machine_work(i) = machine_work(i) + (machineTable{i}(j).end - machineTable{i}(j).start);
        end
    end
end

EG_M_SUM = machineEnergy.work(1: machineNum)' * machine_work + machineEnergy.free(1: machineNum)' * machine_spare;

% [2]. AGV�ܺ�
% ͨ���������ɼ���
EG_AGV = zeros(1, AGVNum);
for i = 1: AGVNum
    for j = 1: size(agvEGRecord{i}, 1)
        if j == 1
            continue;
        end

        if agvEGRecord{i}(j - 1, 2) - agvEGRecord{i}(j, 2) < 0
            continue;
        end

        EG_AGV(i) = EG_AGV(i) + agvEGRecord{i}(j - 1, 2) - agvEGRecord{i}(j, 2);
    end
end

EG_A_SUM = sum(EG_AGV);

FUNC = {[makespan, EG_M_SUM + EG_A_SUM]};

end