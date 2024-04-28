function [machineTable, AGVTable, jobCompleteUnLoad, agvEGRecord, agvChargeNum] = sorting(chrom, jobNum, jobInfo, operaVec, ...
    AGVNum, AGVSpeed, ...
    candidateMachine, distance_matrix, AGVEnergy, AGVEG_MAX, AGVEG_MIN, eChargeSpeed, machineTable, AGVTable)
%% Ⱦɫ��� OS MS AS SS ��Ƭ
operaNum = sum(operaVec);%�洢���й����ܹ�����
OS = chrom(1: operaNum);
MS = chrom(operaNum + 1: 2 * operaNum);
AS = chrom(2 * operaNum + 1 : 3 * operaNum);
SS = chrom(3 * operaNum + 1: 5 * operaNum);

%% ��ѭ��������ѭ��ÿ��Ⱦɫ��
operaRec = zeros(1, jobNum);        % ��¼��ǰ������ÿ��Job�ĵڼ���
curJobTime = zeros(1, jobNum);      % ��¼��ǰJob����һ������Ľ���ʱ��
jobPosition = -1 * ones(1, jobNum); % ��¼��ǰJob����һ������Ļ���
jobCompleteUnLoad = zeros(1, jobNum);   % ��¼ÿ��Job�����͵�ж��վ��ʱ��
agvRealTimeEG = ones(1, AGVNum) * AGVEG_MAX;    % ��¼AGV�ĵ��� ==> ����Լ��
% ������¼�����仯
agvEGRecord = cell(1, AGVNum);
for a = 1: AGVNum
    agvEGRecord{a} = [0, AGVEG_MAX];
end
% ��¼������
agvChargeNum = zeros(1, AGVNum);

%% ѭ��
for i = 1 : sum(operaVec)%��ÿ������
    curJob = OS(i); % ��ǰ������curJob
    operaRec(curJob) = operaRec(curJob) + 1;
    jobOpera = operaRec(curJob);    % ��ǰ�����Ĺ���
    rSIndex = sum(operaVec(1: curJob - 1)) + jobOpera;     % ȷ����ǰjob�ĵ�ǰ�����ӦMS/AS�ϵ�λ��
    machine = candidateMachine{curJob, jobOpera}(MS(rSIndex));  % ��ѡ����
    agv = AS(rSIndex);                                          % ��ѡAGV
    
    % ѡ���ٶȵ�λ
    ssIndex = 2 * (sum(operaVec(1: curJob - 1)) + jobOpera); % �ٶ�����
    freeIndex = ssIndex - 1;
    loadIndex = ssIndex;
    free_speed = AGVSpeed(SS(freeIndex));   % ����SSƬ�ζ�ȡ�������ٶ�
    load_speed = AGVSpeed(SS(loadIndex));   % ����SSƬ�ζ�ȡ�������ٶ�
    freeEGConsume = AGVEnergy.free(SS(freeIndex));  % �����ٶȶ�Ӧ�ĵ�λʱ���ܺ�
    loadEGConsume = AGVEnergy.load(SS(loadIndex));  % �����ٶȶ�Ӧ�ĵ�λʱ���ܺ�


    
    %% AGV������ΪMTSP�������������⣩������˳��ΪOS˳��
    % ���ڲ��ܴ�����ȥ���׮��ÿ��ת�ƣ�����+���أ�ǰ������
    % �������ڳ����ֵ
    % �������AGV

    for ag_idx = 1: AGVNum
        if agvRealTimeEG(ag_idx) <= AGVEG_MIN
            start_m = AGVTable{ag_idx}(end - 1).to_machine;    % AGVĿǰ���ڻ���
            to_m = -2;                                      % ���׮��ж��վ������ -2

            % agv������ж��վ����Ҫǰ��ж��վ���
            if start_m ~= -2
                start_T = AGVTable{ag_idx}(end - 1).end;       % ת�ƿ�ʼʱ��
                transfer_time = distance_matrix.machine_to_unload(start_m) / AGVSpeed(3);    % ת�ƺ�ʱ��ת�Ʋ����ܴ�װ��վ��ʼ
                % ͬʱĬ��ʹ�ÿ���3������쵲λ������Ϊ���벿��û�г�����

                % ǰ��ж��վ���
                insertion.start = start_T;
                insertion.end = insertion.start + transfer_time;
                table_out = table_insert(insertion, AGVTable{ag_idx}, length(AGVTable{ag_idx}), 0, 0, -1, to_m, 2);
                AGVTable{ag_idx} = table_out;

                % ���µ���
                agvRealTimeEG(ag_idx) = agvRealTimeEG(ag_idx) - transfer_time * AGVEnergy.free(3);
                agvEGRecord{ag_idx} = [agvEGRecord{ag_idx}; [insertion.end, agvRealTimeEG(ag_idx)]];
            end

            % ���
            charge_time = (AGVEG_MAX - agvRealTimeEG(ag_idx)) / eChargeSpeed;      % �����������ʱ
            start_T = AGVTable{ag_idx}(end - 1).end;   % ��ʼ���ʱ��
            insertion.start = start_T;
            insertion.end = insertion.start + charge_time;  % �������ʱ��
            table_out = table_insert(insertion, AGVTable{ag_idx}, length(AGVTable{ag_idx}), 0, 0, 0, to_m, 1);
            AGVTable{ag_idx} = table_out;

            % �����󣬵����ָ�
            agvRealTimeEG(ag_idx) = AGVEG_MAX;
            agvEGRecord{ag_idx} = [agvEGRecord{ag_idx}; [insertion.end, agvRealTimeEG(ag_idx)]];
            % ������ +1
            agvChargeNum(ag_idx) = agvChargeNum(ag_idx) + 1;
        end
    end

    agv_complete = curJobTime(curJob);  % ��¼AGV�������γɵ�ʱ��Լ��

    %% ########################### ����ת�� ############################
    len_ = length(AGVTable{agv});
    agv_spare_start_time = AGVTable{agv}(len_).start;
    agv_spare_start_machine = AGVTable{agv}(len_).from_machine;     % ����ת�Ƶ���ʼ����
    agv_spare_dest_machine = jobPosition(curJob);                   % ����ת�Ƶ�Ŀ�����
    
    
    % AGVTable{agv}�������ת��ʱ���
    if agv_spare_dest_machine ~= machine    % ��һ����������һ������ļӹ�������ͬ

        % ����ת��ʱ��
        spare_transfer_time = spare_transfer_time_compute(agv_spare_start_machine, agv_spare_dest_machine, ...
            distance_matrix, free_speed);

        if spare_transfer_time > 1E-6   % spare_transfer_time ~= 0
            insertion.start = agv_spare_start_time;                 % ����ת�Ƶ���ʼʱ��
            insertion.end = insertion.start + spare_transfer_time;  % ����ת�ƵĽ���ʱ��
            table_out = table_insert(insertion, AGVTable{agv}, length(AGVTable{agv}), curJob, jobOpera, -1, ...
                agv_spare_dest_machine, 0);
            AGVTable{agv} = table_out;
            % ���� ����ת�ƺ�AGV�ĵ���
            agvRealTimeEG(agv) = agvRealTimeEG(agv) - spare_transfer_time * freeEGConsume;
            agvEGRecord{agv} = [agvEGRecord{agv}; [insertion.end, agvRealTimeEG(agv)]];
        end
    end

    % 
    %% ########################### ����ת�� ############################
    len_ = length(AGVTable{agv});
    % ����ת�Ƶ���ʼʱ�䣺MAX�������˹��������򣩵��깤ʱ�䡢AGV����ʱ�䣩
    agv_load_start_time = max(curJobTime(curJob), AGVTable{agv}(len_).start);
    agv_load_start_machine = AGVTable{agv}(len_).from_machine;  % ����ת�Ƶ���ʼ����
    agv_load_dest_machine = machine;                            % ����ת�Ƶ�Ŀ�����

    % AGVTable{agv}���븺��ת��ʱ���
    if agv_spare_dest_machine ~= machine    % ��һ����������һ������ļӹ�������ͬ
        
        % ����ת��ʱ��
        load_transfer_time = load_transfer_time_compute(agv_load_start_machine, agv_load_dest_machine, ...
            distance_matrix, load_speed);

        if load_transfer_time > 1E-6   % spare_transfer_time == 0
            insertion.start = agv_load_start_time;                 % ����ת�Ƶ���ʼʱ��
            insertion.end = insertion.start + load_transfer_time;  % ����ת�ƵĽ���ʱ��
            table_out = table_insert(insertion, AGVTable{agv}, length(AGVTable{agv}), curJob, jobOpera, -2, ...
                agv_load_dest_machine, 0);
            AGVTable{agv} = table_out;
            % ���� agv_complete
            agv_complete = insertion.end;
            % ���� ����ת�ƺ�AGV�ĵ���
            agvRealTimeEG(agv) = agvRealTimeEG(agv) - load_transfer_time * loadEGConsume;
            agvEGRecord{agv} = [agvEGRecord{agv}; [insertion.end, agvRealTimeEG(agv)]];
        end
    end

    %% machineTable Ѱ�ղ��빤��ӹ�ʱ����
    for j = 1 : length(machineTable{machine})
        % Ѱ�ҿ���ʱ���
        if isequal(machineTable{machine}(j).job, 0)
            % �ж��ܷ����
            startT = max(machineTable{machine}(j).start, agv_complete);
            endT = startT + jobInfo{curJob}(jobOpera, machine);
            if endT <= machineTable{machine}(j).end           % �ܲ���
                insertion.start=startT;     % ����ʱ������ʼʱ��
                insertion.end=endT;         % ����ʱ���Ľ���ʱ��
                table_out = table_insert(insertion, machineTable{machine}, j, curJob, jobOpera);
                machineTable{machine} = table_out;

                % ����curJobTime��jobPosition 
                curJobTime(curJob) = endT;
                jobPosition(curJob) = machine;
                break;
            end
        end
    end

    %% ����ÿ�����������һ�����򣬰���һ����ʱ��̵�AGV�������䵽 ж��վ
    if jobOpera == operaVec(curJob)
        arrival_time = [];
        for a = 1 : AGVNum
            earlierest_start_time = AGVTable{a}(length(AGVTable{a})).start;
            earlierest_start_machine = AGVTable{a}(length(AGVTable{a})).from_machine;
            % �����ӹ����������ж��վ���˲����ٶ�û�б��룬Ĭ�ϲ�������
            return_free_speed = AGVSpeed(3);
            earlierest_transfer_time = spare_transfer_time_compute(earlierest_start_machine, machine, distance_matrix, return_free_speed);
            earlierest_arrival_time = earlierest_start_time + earlierest_transfer_time;
            arrival_time = [arrival_time, earlierest_arrival_time];
        end

        % ѡ����������뿪��
        leave_time = max([ones(1, AGVNum) * curJobTime(curJob); arrival_time], [], 1);
        agv_candidate = find(leave_time == min(leave_time));
        if length(agv_candidate) > 1
            new_arrival_time = arrival_time(agv_candidate);
            last_arrival_index = find(new_arrival_time == max(new_arrival_time));
            return_agv = agv_candidate(last_arrival_index(1));
        else
            return_agv = agv_candidate;
        end

        % ����ǰ�����˹���
        if AGVTable{return_agv}(length(AGVTable{return_agv})).from_machine ~= machine
            insertion.start = AGVTable{return_agv}(length(AGVTable{return_agv})).start;
            insertion.end = arrival_time(return_agv); 
            table_out = table_insert(insertion, AGVTable{return_agv}, length(AGVTable{return_agv}), ...
                curJob, -1, -1, machine, 0);
            AGVTable{return_agv} = table_out;
            % �����������µ��� return_agv
            agvRealTimeEG(return_agv) = agvRealTimeEG(return_agv) - (insertion.end - insertion.start) * AGVEnergy.free(3);
            agvEGRecord{return_agv} = [agvEGRecord{return_agv}; [insertion.end, agvRealTimeEG(return_agv)]];
        end

        % returnAGV ���ù��������� ж��վ
        return_load_start_time = max(arrival_time(return_agv), curJobTime(curJob));
        % �����ӹ����������ж��վ���˲����ٶ�û�б��룬Ĭ�ϲ�������
        return_load_speed = AGVSpeed(3);
        return_load_transfer_time = distance_matrix.machine_to_unload(machine) / return_load_speed;
        insertion.start = return_load_start_time;
        insertion.end = insertion.start + return_load_transfer_time;
        table_out = table_insert(insertion, AGVTable{return_agv}, length(AGVTable{return_agv}), ...
                curJob, -1, -2, -2, 0);
        AGVTable{return_agv} = table_out;
        % �����������µ��� return_agv
        agvRealTimeEG(return_agv) = agvRealTimeEG(return_agv) - (insertion.end - insertion.start) * AGVEnergy.load(3);
        agvEGRecord{return_agv} = [agvEGRecord{return_agv}; [insertion.end, agvRealTimeEG(return_agv)]];

        % ���� jobCompleteUnLoad
        jobCompleteUnLoad(curJob) = insertion.end;

    end
    
end
end