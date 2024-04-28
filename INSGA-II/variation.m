function off_spring = variation(p_cross, p_mutation, parent_chromosome, jobNum, operaVec, AGVNum, AGVSpeed, candidateMachine)

N_size = size(parent_chromosome, 1);
operaNum = sum(operaVec);       % �ܹ�����
n_var = 5 * operaNum;           % OS��MS AS SS ��������
os_len = sum(operaVec);         % OS��������
rs_len = 4 * sum(operaVec);     % RS��������  RS = MS + AS + SS
off_spring=[];

% UP LOW �������   UP��RS�Ͻ�    LOW��RS�½�
UP = [];
for i = 1: jobNum
    for j = 1: operaVec(i)
        % MS���Ͻ�
        UP = [UP length(candidateMachine{i, j})];
    end
end
UP = [UP AGVNum * ones(1, operaNum) length(AGVSpeed) * ones(1, 2 * operaNum)];

%% ���������ѭ��
for i = 1: N_size
    parent_1 = round(N_size * rand);
    if parent_1 < 1
        parent_1 = 1;
    end
    child = parent_chromosome(parent_1, 1: n_var);
    
    %% ����
    % OS ==> IPOX����
    % RS ==> MPX����
    if rand < p_cross
        child = [];
        parent_2 = round(N_size * rand);
        if parent_2 < 1
            parent_2 = 1;
        end
        while isequal(parent_chromosome(parent_1, 1: n_var), parent_chromosome(parent_2, 1: n_var))
            parent_2 = round(N_size * rand);
            if parent_2 < 1
                parent_2 = 1;
            end
        end

        parent_1_os = parent_chromosome(parent_1, 1: os_len);           % P1 [OS]
        parent_1_rs = parent_chromosome(parent_1, os_len + 1: n_var);   % P1 [MS AS SS]
        parent_2_os = parent_chromosome(parent_2, 1: os_len);           % P2 [OS]
        parent_2_rs = parent_chromosome(parent_2, os_len + 1: n_var);   % P2 [MS AS SS]
        
        %% IPOX����
        % OS ����
        random_r = randperm(jobNum, 1);
        job_set = sort(randperm(jobNum, random_r));
        sequence_1_index = ismember(parent_1_os, job_set);
        sequence_2_index = ismember(parent_2_os, job_set);
        child_1_os = parent_1_os .* ~ sequence_1_index;
        child_1_os(child_1_os == 0) = parent_2_os(sequence_2_index);
        child_2_os = parent_2_os .* ~ sequence_2_index;
        child_2_os(child_2_os == 0) = parent_1_os(sequence_1_index);
        
        % RS��[MS AS SS]���֣�MPX����
        random_r = randperm(rs_len, 1);
        random_r = sort(randperm(rs_len, random_r));
        child_2_rs = parent_1_rs;
        child_1_rs = parent_2_rs;
        for k = 1: length(random_r)
            child_1_rs(random_r(k)) = parent_1_rs(random_r(k));
            child_2_rs(random_r(k)) = parent_2_rs(random_r(k));
        end

        % OS��RS���
        child(1, :) = [child_1_os, child_2_rs];
        child(2, :) = [child_2_os, child_1_rs];
    end
    
    %% �Ƿ���б���
    for child_num = 1: size(child, 1)
        if rand < p_mutation
            % OS ==> ��������
            posi1 = randperm(os_len, 1);
            posi2 = randperm(os_len, 1);
            while child(child_num, posi1) == child(child_num, posi2)
                posi2 = randperm(os_len, 1);
            end
            position = [posi1 posi2];
            os_selection = child(child_num, position);
            os_cache = os_selection(1);
            os_selection(1) = os_selection(2);
            os_selection(2) = os_cache;
            child(child_num, position) = os_selection;

            % RS ==> ������
            random_r = max(1, randperm(round(0.05 * rs_len), 1));
            random_r = randperm(rs_len, random_r);
            for k = 1: length(random_r)
                child(child_num, os_len + random_r(k)) = randperm(UP(random_r(k)), 1);
            end
            
        end
    end

    for j = 1: size(child, 1)
        off_spring(size(off_spring, 1) + 1, :) = child(j, :);
    end
end
end