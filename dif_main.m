clear
clc
close all

numrun = 1;
for i = 1:numrun
    % clear
    clc
    close all
    p_cross = 0.8;              % 交叉概率
    p_mutation = 0.2;           % 变异概率
    min_pr = 0.05;              % 反向生成比例下限
    max_pr = 0.3;              % 反向生成比例上限
    epsilon = 0.8;
    alpha = 0.1;
    gamma = 0.9;


    pth = 'fjsp\Brandimarte_Data\Mk02.fjs';%2,4,6,8
    benchmarkRead(pth);
    load data.mat

    distance_from_xy(machineNum);
    distance_matrix_excel = xlsread('机器数据.xlsx', '装卸站到机器距离');
    distance_matrix.load_to_machine = distance_matrix_excel(1, :);
    distance_matrix.load_to_machine = distance_matrix.load_to_machine(1: machineNum);
    distance_matrix.machine_to_unload = distance_matrix_excel(2, :);
    distance_matrix.machine_to_unload = distance_matrix.machine_to_unload(1: machineNum);
    distance_matrix.machine_to_machine = xlsread('机器数据.xlsx', '机器到机器距离');
    distance_matrix.machine_to_machine = distance_matrix.machine_to_machine(1: machineNum, 1: machineNum);
    distance_matrix.load_to_unload = xlsread('机器数据.xlsx', '装载站到卸载站距离');
    machineEnergy.work = xlsread('机器数据.xlsx', '机器加工能耗');
    machineEnergy.free = xlsread('机器数据.xlsx', '机器空载能耗');

    AGVNum = 3;
    AGVSpeed = [0.5,0.75,1.0];
    % AGV能耗
    AGVEnergy_excel = xlsread('AGV数据.xlsx', 'AGV能耗');
    AGVEnergy.free = AGVEnergy_excel(1, :);
    AGVEnergy.load = AGVEnergy_excel(2, :);
    clear distance_matrix_excel AGVEnergy_excel

    AGVEG_MIN = 19.2;
    AGVEG_MAX = 100;

    distance_MAX = max([max(distance_matrix.machine_to_machine) ...
        distance_matrix.load_to_machine ...
        distance_matrix.machine_to_unload...
        distance_matrix.load_to_unload]);
    check_MIN = 1.0 * distance_MAX / AGVSpeed(end) * (AGVEnergy.free(end) + AGVEnergy.load(end));
    disp(['power > ' num2str(check_MIN)])
    if check_MIN > AGVEG_MIN
        error('power error')
    end
    eChargeSpeed = 20;
    pop = 10;
    max_gen = 10;
    speedNum = length(AGVSpeed); 
    operaNum = sum(operaNumVec); 
    chrom = init(pop, jobNum, operaNumVec, candidateMachine, AGVNum, speedNum);
    cd('INSGA-II\')
    INSGA_II_0_Result = INSGA_II(p_cross,p_mutation,min_pr,max_pr,epsilon,alpha,gamma,pop,chrom,max_gen,jobNum, jobInfo, operaNumVec, machineNum, AGVNum, AGVSpeed, candidateMachine, distance_matrix, machineEnergy, ...
        AGVEnergy, AGVEG_MAX, AGVEG_MIN, eChargeSpeed, 'VNS+i-elitism');
    cd('..\')
end
beep;