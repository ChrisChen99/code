  function NSGA2_Result = INSGA_II(p_cross,p_mutation,min_pr,max_pr,epsilon,alpha,gamma,pop,chrom,max_gen,jobNum, jobInfo, operaVec, machineNum, AGVNum, AGVSpeed, candidateMachine, distance_matrix, machineEnergy, ...
    AGVEnergy, AGVEG_MAX, AGVEG_MIN, eChargeSpeed, strategy)
tstart = tic;               % ��ʱ��
disp(['RUNNING --------> NSGA-II ' strategy ' <-------- RUNNING'])
disp(['job��' num2str(jobNum), ' machine ', num2str(machineNum), ' AGV ' num2str(AGVNum)]);

% ģ����Ϣ
speedNum = length(AGVSpeed);    % AGV�ٶȵ�λ��Ŀ
operaNum = sum(operaVec);       % ��������

%% NSGA-II ����
dim = 5 * operaNum;         
obj_num = 2;                
pop = 150;                  
pool = round(pop / 2);     
tour = 2;                  

gen = 0;
Qtable=ones(4,8);
currt_State = 4;
next_State = 1;
last_HV = 0;
last_SP = 0;
reward = 0;
pointer = 0;
%% ��Ⱥ��ʼ��
chrom = init(pop, jobNum, operaVec, candidateMachine, AGVNum, speedNum);
for i = 1 : pop
    func = fitness(chrom(i, :), jobNum, jobInfo, operaVec, machineNum, AGVNum, AGVSpeed, candidateMachine, ...
        distance_matrix, machineEnergy, AGVEnergy, AGVEG_MAX, AGVEG_MIN, eChargeSpeed);
    chrom(i, dim + 1: dim + obj_num) = func{1};
end
pop_history = cell(1, max_gen);
%% ������ѭ��
while gen < max_gen
    parent_ = tournament_selection(chrom, pool, tour);
    offspring_ = variation(p_cross, p_mutation, parent_, jobNum, operaVec, AGVNum, AGVSpeed, candidateMachine);
    for i = 1: size(offspring_, 1)
        func = fitness(offspring_(i, :), jobNum, jobInfo, operaVec, machineNum, AGVNum, AGVSpeed, candidateMachine, ...
            distance_matrix, machineEnergy, AGVEnergy, AGVEG_MAX, AGVEG_MIN, eChargeSpeed);
        offspring_(i, dim + 1: dim + obj_num) = func{1};
    end
    intermediate_ = [chrom(:, 1: dim + obj_num); offspring_(:, 1: dim + obj_num)];
    intermediate_ = non_domination(intermediate_, dim, obj_num);%�Ӹ���2N
    elitism_pop = intermediate_(1: pop, :);%ѡ��ǰN�����뾫Ӣ��Ⱥ

    %% ��1������ѧϰ
    if contains(strategy, 'i-elitism')
            reverse_percent = min_pr + (max_pr - min_pr) * (max_gen - gen) / (max_gen);
        elitism_pop = improved_elitism(elitism_pop, reverse_percent, jobNum, jobInfo, operaVec, machineNum, AGVNum, AGVSpeed, candidateMachine, ...
            distance_matrix, machineEnergy, AGVEnergy, AGVEG_MAX, AGVEG_MIN, eChargeSpeed, obj_num);
        %end
    end
    %% ��2���ֲ�����
    chrom = elitism_pop;
    vns_chrom = [];
    if contains(strategy, 'VNS')
        if gen >=5 && isequal(curve_min(:, gen), curve_min(:, gen-1)) &&...
                isequal(curve_min(:, gen), curve_min(:, gen-2))&&...
                isequal(curve_min(:, gen), curve_min(:, gen-3))
            pointer = 1;
        end
        if pointer == 1
            option_ = {'O1+M1+A1', 'O1+M1+A2', 'O1+M2+A1', 'O1+M2+A2', ...
                'O2+M1+A1', 'O2+M1+A2', 'O2+M2+A1', 'O2+M2+A2'};
            %ѡ��ǰi���ֲ���������
            if rand>epsilon||all(Qtable(currt_State,:)==1)%%e-greedy����
                action = randperm(numel(option_), 1);%%����һ�������������i������
            else
                [~,action] = max(Qtable(currt_State,:));
            end
            for i = 1: size(chrom, 1)
                vns_ = [];
                vns_ = VNS(chrom(i, 1: dim + obj_num), jobNum, jobInfo, operaVec, machineNum, AGVNum, AGVSpeed, candidateMachine, ...
                    distance_matrix, machineEnergy, AGVEnergy, AGVEG_MAX, AGVEG_MIN, eChargeSpeed, obj_num, option_{action});
                vns_chrom = [vns_chrom; vns_];
            end

            intermediate_ = [chrom(:, 1: dim + obj_num); vns_chrom];
            intermediate_ = non_domination(intermediate_, dim, obj_num);
            elitism_pop = intermediate_(1: pop, :);
            chrom = elitism_pop;

            first_front_individual = chrom(chrom(:, dim + obj_num + 1) == 1, :);
            first_valueset  = first_front_individual(:, dim + 1:dim + obj_num);
            ref_point = [160.0; 2100.0]; 
            ref_point = double(ref_point);
            
            cd('HV\')
            currt_HV = test_lebesgue_measure(first_valueset,ref_point);
            cd('..\')
            % Spacingָ��
            cd('Spacing\')
            currt_SP =  Spacing(first_valueset);
            cd('..\')

            delt_HV=currt_HV-last_HV;
            delt_SP=currt_SP-last_SP;

            if delt_HV>0&& delt_SP<0
                reward = 2;
            elseif delt_HV>0&& delt_SP>=0
                reward = 0;
            elseif delt_HV<=0&& delt_SP<0
                reward = 0;
            elseif delt_HV<=0&& delt_SP>=0
                reward = -1;
            end

            next_State = getState(delt_HV,delt_SP);
            maxreward = max(Qtable(next_State,:));

            Qtarget=reward+gamma*maxreward;
            Qtable(currt_State,action)=Qtable(currt_State,action)+...
                alpha*(Qtarget-Qtable(currt_State,action));
            %����Q����

            last_HV = currt_HV;
            last_SP = currt_SP;
            currt_State = next_State;
            reward = 0;  
        end
    end


    % ����ǰ������Ⱥ���浽 cell ������
    pop_history{gen + 1} = chrom;
    gen = gen + 1;

    %% Disp
    fprintf('GEN: %d  MIN Cmax: %.1f  MIN Energy:%.2f\n', gen, min(chrom(:, dim + 1)), min(chrom(:, dim + 2)));
    for k = 1: obj_num
        curve_min(k, gen) = min(chrom(:, dim + k));
        curve_avg(k, gen) = mean(chrom(:, dim + k));
    end
end
% ����Ⱥ��ʷ��Ϣ���浽����ṹ��
NSGA2_Result.pop_history = pop_history;
RunTime = toc(tstart);
disp(['����ʱ�䣺' num2str(RunTime)]);

%%  ����֧����¼
chrom = chrom(chrom(:, dim + obj_num + 1) == 1, :);        % ֧��ȼ�Ϊ1�ĸ���
obj_matrix = chrom(:, dim + 1: dim + obj_num);
[obj_matrix, uni_index] = unique(obj_matrix, 'rows');
chrom = chrom(uni_index, :);

%% ����
NSGA2_Result.RunTime = RunTime;
NSGA2_Result.chrom = chrom;
NSGA2_Result.obj_matrix = obj_matrix;
NSGA2_Result.curve.min = curve_min;
NSGA2_Result.curve.avg = curve_avg;
for idx = 1: size(chrom, 1)
    [~, NSGA2_Result.machineTable{idx}, NSGA2_Result.AGVTable{idx}, ~, NSGA2_Result.EG_M_SUM{idx}, ...
        NSGA2_Result.EG_A_SUM{idx}, NSGA2_Result.agvEGRecord{idx}, NSGA2_Result.agvChargeNum{idx}]...
        = fitness(chrom(idx, 1: dim), jobNum, jobInfo, operaVec, machineNum, AGVNum, AGVSpeed, candidateMachine, ...
        distance_matrix, machineEnergy, AGVEnergy, AGVEG_MAX, AGVEG_MIN, eChargeSpeed);
end