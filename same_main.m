clc
clear
close all
% 指定用于保存.fig图形的文件夹
figureSaveFolder = 'figures';

if ~exist(figureSaveFolder, 'dir')
    mkdir(figureSaveFolder);
end
numrun = 2;
% param_settings = [
% %     0.7, 0.1, 0.1,0.4,0.7,0.1,0.7;%3
% %     0.7, 0.2, 0.2,0.5,0.8,0.2,0.8;%1(最优)
% %     0.7, 0.3, 0.3,0.6,0.9,0.3,0.9;%2
% %     0.8, 0.1, 0.2,0.6,0.7,0.2,0.9;
% %     0.8, 0.2, 0.3,0.4,0.8,0.3,0.7;
% %     0.8, 0.3, 0.1,0.5,0.9,0.1,0.8;
%     0.9, 0.1, 0.3,0.5,0.7,0.3,0.8;
%     0.9, 0.2, 0.1,0.6,0.8,0.1,0.9;
%     0.9, 0.3, 0.2,0.4,0.9,0.2,0.7;%%%%
%     ];
%
% resultCell_I = cell(size(param_settings, 1), 1);
% %resultCell_NSGA = cell(size(param_settings, 1), 1);
% for j = 1:size(param_settings, 1)
%     p_cross = param_settings(j, 1);             % 交叉概率
%     p_mutation = param_settings(j, 2);          % 变异概率
%     min_pr = param_settings(j, 3);     % 反向生成比例
%     max_pr = param_settings(j, 4);     % 反向生成比例
%     epsilon = param_settings(j, 5);
%     alpha = param_settings(j, 6);
%     gamma = param_settings(j, 7);
%
%
%     numresult_I = [];
%     %numresult_NSGA = [];
for i = 1:numrun
    clc
    close all
    disp(['当前是第 ' num2str(i) ' 次独立运行']);
    %
    %         disp(['第' num2str(j) '组参数的第' num2str(i) '次独立运行']);

    p_cross = 0.8;              % 交叉概率
    p_mutation = 0.2;           % 变异概率
    min_pr =0.1;              % 反向生成比例下限
    max_pr = 0.3;              % 反向生成比例上限
    epsilon = 0.8;
    alpha = 0.1;
    gamma = 0.8;

    %% 载入标准算例的时间数据
    pth = 'fjsp\Brandimarte_Data\Mk02.fjs';%2,4,8
    %pth = 'fjsp\Dauzere_Data\18a.fjs';%4,8
    %pth = 'fjsp\Barnes\mt10x.fjs';%4,8
    benchmarkRead(pth);
    load data.mat

    %% 与机器/装卸载站相关数据
    %
    % 使用坐标（x, y）生成距离数据
    distance_from_xy(machineNum);
    distance_matrix_excel = xlsread('机器数据.xlsx', '装卸站到机器距离');
    % 装载站到每台机器的·距离
    distance_matrix.load_to_machine = distance_matrix_excel(1, :);
    distance_matrix.load_to_machine = distance_matrix.load_to_machine(1: machineNum);
    % 卸载站到每台机器的距离
    distance_matrix.machine_to_unload = distance_matrix_excel(2, :);
    distance_matrix.machine_to_unload = distance_matrix.machine_to_unload(1: machineNum);
    % 每台机器间的距离
    distance_matrix.machine_to_machine = xlsread('机器数据.xlsx', '机器到机器距离');
    distance_matrix.machine_to_machine = distance_matrix.machine_to_machine(1: machineNum, 1: machineNum);
    % 装 卸载站之间的距离
    distance_matrix.load_to_unload = xlsread('机器数据.xlsx', '装载站到卸载站距离');
    % 能耗
    machineEnergy.work = xlsread('机器数据.xlsx', '机器加工能耗');
    machineEnergy.free = xlsread('机器数据.xlsx', '机器空载能耗');

    %% AGV相关数据
    AGVNum = 3;
    %AGVSpeed = [0.5,0.75,1];
    AGVSpeed = [0.5,0.75,1];
    % AGV能耗
    AGVEnergy_excel = xlsread('AGV数据.xlsx', 'AGV能耗');
    AGVEnergy.free = AGVEnergy_excel(1, :);
    AGVEnergy.load = AGVEnergy_excel(2, :);

    clear distance_matrix_excel AGVEnergy_excel

    %% 充电桩数据
    % [1]. 充电桩放置在 卸载站 ==> 充电桩位置信息与卸载站相同
    % AGV额定电量 Kw*h
    AGVEG_MAX = 100;
    % AGV充电阈值 Kw*h ↙↙↙↙↙
    % mk01 20kW*h
    % mk02 19.2kW*h
    % mk03 kW*h
    % mk04 24kW*h
    % mk08 24kW*h
    AGVEG_MIN = 19.2;
    % 验证充电阈值是否合理
    distance_MAX = max([max(distance_matrix.machine_to_machine) ...
        distance_matrix.load_to_machine ...
        distance_matrix.machine_to_unload...
        distance_matrix.load_to_unload]);
    % 电量最低阈值
    % 1.1表示多增加 10% 电量
    check_MIN = 1.0 * distance_MAX / AGVSpeed(end) * (AGVEnergy.free(end) + AGVEnergy.load(end));
    disp(['电量阈值 > ' num2str(check_MIN)])
    if check_MIN > AGVEG_MIN
        error('AGV充电阈值 ==> 设置错误')
    end
    % 充电速度 Kw
    eChargeSpeed = 20;

    %% NSAG-II 算法
    % jobNum            工件数量
    % jobInfo           工件加工时间信息
    % operaNumVec       每个工件的工序数量
    % machineNum        机器数量
    % AGVNum            AGV数量
    % AGVSpeed          AGV速度（不同挡位）
    % candidateMachine  每个工件的候选加工机器
    % distance_matrix   距离矩阵
    % machineEnergy     机器能耗
    % AGVEnergy         AGV能耗（不同挡位）
    % AGVEG_MAX         AGV额定电量
    % AGVEG_MIN         AGV充电阈值（低于此电量，需进行充电）
    % eChargeSpeed      AGV充电速率
    pop = 100;
    max_gen = 200;
    speedNum = length(AGVSpeed);    % AGV速度位数目
    operaNum = sum(operaNumVec);       % 工序总数
    %%种群初始化
    chrom = init(pop, jobNum, operaNumVec, candidateMachine, AGVNum, speedNum);

    % 算法对比
    %%%【1】NSGA
    % cd('NSGA-II\')
    % NSGA_II_Result = NSGA2(p_cross,p_mutation,pop,chrom,max_gen,jobNum, jobInfo, operaNumVec, machineNum, AGVNum, AGVSpeed, candidateMachine, distance_matrix, machineEnergy, ...
    %     AGVEnergy, AGVEG_MAX, AGVEG_MIN, eChargeSpeed);
    % cd('..\')
    %
    % % 【3】反向
    % cd('INSGA-II\')
    % INSGA_II_1_Result = INSGA_II(p_cross,p_mutation,min_pr,max_pr,epsilon,alpha,gamma,pop,chrom,max_gen,jobNum, jobInfo, operaNumVec, machineNum, AGVNum, AGVSpeed, candidateMachine, distance_matrix, machineEnergy, ...
    %     AGVEnergy, AGVEG_MAX, AGVEG_MIN, eChargeSpeed, 'i-elitism');
    % cd('..\')

    % % %【1】i-elitism + 不加Q-Learning，VNS
    % cd('INSGA-II\')
    % INSGA_II_3_Result = INSGA_II(p_cross,p_mutation,min_pr,max_pr,epsilon,alpha,gamma,pop,chrom,max_gen,jobNum, jobInfo, operaNumVec, machineNum, AGVNum, AGVSpeed, candidateMachine, distance_matrix, machineEnergy, ...
    %     AGVEnergy, AGVEG_MAX, AGVEG_MIN, eChargeSpeed, 'NOQ_VNS+i-elitism');
    % cd('..\')
    %
    %【1】i-elitism + VNS
    cd('INSGA-II\')
    INSGA_II_0_Result = INSGA_II(p_cross,p_mutation,min_pr,max_pr,epsilon,alpha,gamma,pop,chrom,max_gen,jobNum, jobInfo, operaNumVec, machineNum, AGVNum, AGVSpeed, candidateMachine, distance_matrix, machineEnergy, ...
        AGVEnergy, AGVEG_MAX, AGVEG_MIN, eChargeSpeed, 'VNS+i-elitism');
    cd('..\')
    %
    % %%【2】VNS
    % cd('INSGA-II\')
    % INSGA_II_2_Result = INSGA_II(p_cross,p_mutation,min_pr,max_pr,epsilon,alpha,gamma,pop,chrom,max_gen,jobNum, jobInfo, operaNumVec, machineNum, AGVNum, AGVSpeed, candidateMachine, distance_matrix, machineEnergy, ...
    %     AGVEnergy, AGVEG_MAX, AGVEG_MIN, eChargeSpeed, 'VNS');
    % cd('..\')


    %
    % figure(2);
    % plot(NSGA_II_Result.curve.min(1, :), 'Color', [0, 0.447, 0.741], 'LineStyle', '-', 'LineWidth', 1.5) % 蓝色
    % hold on
    % plot(INSGA_II_2_Result.curve.min(1, :), 'Color', [0.85, 0.325, 0.098], 'LineStyle', '-', 'LineWidth', 1.5) % 红色
    % hold on
    % grid on; box on;
    % legend({'NSGA-II', 'I-NSGA-II-ML'}, 'FontName', 'Times New Roman')
    % xlabel('迭代次数', 'FontName', 'Times New Roman')
    % ylabel('最大完工时间 \ith', 'FontName', 'Times New Roman')
    % set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)
    % set(gcf, 'Position', [100, 100, 800, 600])  % 调整图的大小


    % % 创建一个图形窗口
    % figure(1);
    % % 曲线数据
    % x1 = NSGA_II_Result.curve.min(1, :);
    % x2 = INSGA_II_2_Result.curve.min(1, :);
    % % 定义颜色
    % colors = {[0, 0.447, 0.741], [0.85, 0.325, 0.098]}; % 蓝色和红色
    % % 绘制两条曲线
    % plot(x1, 'Color', colors{1}, 'LineStyle', '-', 'LineWidth', 1.5);
    % hold on;
    % plot(x2, 'Color', colors{2}, 'LineStyle', '-', 'LineWidth', 1.5);
    % % 标记每10代，从第5代开始
    % indices = 5:10:numel(x1);
    % for j = 1:length(indices)
    %     plot(indices(j), x1(indices(j)), 's', 'MarkerSize', 10, 'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', 'k');
    %     plot(indices(j), x2(indices(j)), 's', 'MarkerSize', 10, 'MarkerFaceColor', colors{2}, 'MarkerEdgeColor', 'k');
    % end
    % grid on;
    % box on;
    % legend({'','','NSGA-II', 'I-NSGA-II-ML'}, 'FontName', 'Times New Roman')
    % % 创建图例
    % xlabel('迭代次数', 'FontName', 'Times New Roman');
    % ylabel('最大完工时间 \ith', 'FontName', 'Times New Roman');
    % set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
    % set(gcf, 'Position', [100, 100, 800, 600]);
    % set(gca,'Box','off');
    %
    % % %% 迭代曲线
    % %
    % % %         fig1 = figure('Name', ['时间-第' num2str(j) '的参数第' num2str(i) '组参数运行']);
    % % %         % 在这里添加绘图或数据操作
    % % figure(1);
    % % plot(NSGA_II_Result.curve.min(1, :), 'blue', 'LineStyle', '--','LineWidth', 1.3)
    % % hold on
    % % plot(INSGA_II_0_Result.curve.min(1, :), 'red', 'LineStyle', '-','LineWidth', 1.3)
    % % hold on
    % % %
    % % %         %测试
    % % %         legend({'NSGA-II', 'I-NSGA-II-ML'}, 'FontName', 'Times New Roman')
    % %
    % % plot(INSGA_II_1_Result.curve.min(1, :), 'magenta','LineStyle', ':', 'LineWidth', 1.3)
    % % hold on
    % % plot(INSGA_II_2_Result.curve.min(1, :), 'black', 'LineStyle', '-.','LineWidth', 1.3)
    % % grid on; box on;
    % % legend({'NSGA-II', 'I-NSGA-II-ML', 'I-NSGA-II-OBL', 'I-NSGA-II-QL'}, 'FontName', 'Times New Roman')
    % % xlabel('迭代次数')
    % % ylabel('最大完工时间 \ith')
    % % title('迭代曲线')
    % %
    % % %
    % % %fig2 = figure('Name', ['能耗-第' num2str(j) '参数的第' num2str(i) '组参数运行']);
    % % % 在这里添加绘图或数据操作
    % % figure(2);
    % % plot(NSGA_II_Result.curve.min(2, :), 'blue', 'LineStyle', '--','LineWidth', 1.3)
    % % hold on
    % % plot(INSGA_II_0_Result.curve.min(2, :), 'red','LineStyle', '-', 'LineWidth', 1.3)
    % % hold on
    % % %%测试
    % % %legend({'NSGA-II', 'I-NSGA-II-ML'}, 'FontName', 'Times New Roman')
    % %
    % % plot(INSGA_II_1_Result.curve.min(2, :), 'magenta', 'LineStyle', ':','LineWidth', 1.3)
    % % hold on
    % % plot(INSGA_II_2_Result.curve.min(2, :), 'black', 'LineStyle', '-.','LineWidth', 1.3)
    % % grid on; box on;
    % % legend({'NSGA-II', 'I-NSGA-II-ML', 'I-NSGA-II-OBL', 'I-NSGA-II-QL'}, 'FontName', 'Times New Roman')
    % % xlabel('迭代次数')
    % % ylabel('总能耗 \itKw')
    % % title('迭代曲线')
    % %
    % %
    % %
    % % figure(3);
    % % plot(NSGA_II_Result.curve.min(1, :), 'blue', 'LineStyle', '--','LineWidth', 1.3)
    % % hold on
    % % plot(INSGA_II_0_Result.curve.min(1, :), 'red','LineStyle', '-', 'LineWidth', 1.3)
    % % hold on
    % % grid on; box on;
    % % legend({'NSGA-II', 'I-NSGA-II-ML'}, 'FontName', 'Times New Roman')
    % % xlabel('迭代次数')
    % % ylabel('最大完工时间 \itKw')
    % % title('迭代曲线')
    % %
    % % figure(4);
    % % plot(NSGA_II_Result.curve.min(2, :), 'blue', 'LineStyle', '--','LineWidth', 1.3)
    % % hold on
    % % plot(INSGA_II_0_Result.curve.min(2, :), 'red','LineStyle', '-', 'LineWidth', 1.3)
    % % hold on
    % % grid on; box on;
    % % legend({'NSGA-II', 'I-NSGA-II-ML'}, 'FontName', 'Times New Roman')
    % % xlabel('迭代次数')
    % % ylabel('总能耗 \itKw')
    % % title('迭代曲线')
    % %
    % % figure(5);
    % % plot(NSGA_II_Result.curve.min(1, :), 'blue', 'LineStyle', '--','LineWidth', 1.3)
    % % hold on
    % % plot(INSGA_II_1_Result.curve.min(1, :), 'magenta','LineStyle', ':', 'LineWidth', 1.3)
    % % hold on
    % % grid on; box on;
    % % legend({'NSGA-II', 'I-NSGA-II-OBL'}, 'FontName', 'Times New Roman')
    % % xlabel('迭代次数')
    % % ylabel('最大完工时间 \itKw')
    % % title('迭代曲线')
    %
    % % figure(5);
    % % plot(NSGA_II_Result.curve.min(1, :), 'b--', 'LineWidth', 1.3)
    % % hold on
    % % plot(INSGA_II_1_Result.curve.min(1, :), 'm:', 'LineWidth', 1.3)
    % % grid on; box on;
    % % legend({'NSGA-II', 'I-NSGA-II-OBL'}, 'FontName', 'Times New Roman', 'FontSize', 12)
    % % xlabel('迭代次数', 'FontName',  'FontSize', 12)
    % % ylabel('最大完工时间 \itKw', 'FontName', 'FontSize', 12)
    % % title('迭代曲线', 'FontName', 'FontSize', 14)
    % % set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)
    % % set(gcf, 'Position', [100, 100, 800, 600])  % 调整图的大小
    % %
    % % figure(5);
    % % % 使用不同的颜色映射
    % % cmap = colormap(parula(2)); % 这里使用了Parula调色板，其中2表示两种颜色
    % % colormap(cmap);
    % % plot(NSGA_II_Result.curve.min(1, :), 'Color', cmap(2, :), 'LineWidth', 1.3) % 指定颜色为cmap(1, :)
    % % hold on
    % % plot(INSGA_II_1_Result.curve.min(1, :), 'Color', cmap(1, :), 'LineWidth', 1.3) % 指定颜色为cmap(2, :)
    % % % 添加数据点标记
    % % scatter(1:length(NSGA_II_Result.curve.min(1, :)), NSGA_II_Result.curve.min(1, :), 40, cmap(2, :), 'v', 'filled')
    % % scatter(1:length(INSGA_II_1_Result.curve.min(1, :)), INSGA_II_1_Result.curve.min(1, :), 40, cmap(1, :), 'v', 'filled')
    % % grid on; box on;
    % % legend({'NSGA-II', 'I-NSGA-II-OBL'}, 'FontName', 'Times New Roman', 'FontSize', 12)
    % % xlabel('迭代次数', 'FontName', '宋体', 'FontSize', 12)
    % % ylabel('最大完工时间 \ith', 'FontName', '宋体', 'FontSize', 12)
    % % title('迭代曲线', 'FontName', '宋体', 'FontSize', 14)
    % % set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)
    % % set(gcf, 'Position', [100, 100, 800, 600])  % 调整图的大小
    % %
    % % figure(6);
    % % % 使用不同的颜色映射
    % % cmap = colormap(parula(2)); % 这里使用了Parula调色板，其中2表示两种颜色
    % % colormap(cmap);
    % % plot(NSGA_II_Result.curve.min(2, :), 'Color', cmap(2, :), 'LineWidth', 1.3) % 指定颜色为cmap(1, :)
    % % hold on
    % % plot(INSGA_II_1_Result.curve.min(2, :), 'Color', cmap(1, :), 'LineWidth', 1.3) % 指定颜色为cmap(2, :)
    % % % 添加数据点标记
    % % scatter(1:length(NSGA_II_Result.curve.min(2, :)), NSGA_II_Result.curve.min(2, :), 40, cmap(2, :), 'v', 'filled')
    % % scatter(1:length(INSGA_II_1_Result.curve.min(2, :)), INSGA_II_1_Result.curve.min(2, :), 40, cmap(1, :), 'v', 'filled')
    % % grid on; box on;
    % % legend({'NSGA-II', 'I-NSGA-II-OBL'}, 'FontName', 'Times New Roman', 'FontSize', 12)
    % % xlabel('迭代次数', 'FontName', '宋体', 'FontSize', 12)
    % % ylabel('总能耗 \itKw', 'FontName', '宋体', 'FontSize', 12)
    % % title('迭代曲线', 'FontName', '宋体', 'FontSize', 14)
    % % set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)
    % % set(gcf, 'Position', [100, 100, 800, 600])  % 调整图的大小
    % % figure(6);
    % % plot(NSGA_II_Result.curve.min(2, :), 'blue', 'LineStyle', '--','LineWidth', 1.3)
    % % hold on
    % % plot(INSGA_II_1_Result.curve.min(2, :), 'magenta','LineStyle', ':', 'LineWidth', 1.3)
    % % hold on
    % % grid on; box on;
    % % legend({'NSGA-II', 'I-NSGA-II-OBL'}, 'FontName', 'Times New Roman')
    % % xlabel('迭代次数')
    % % ylabel('总能耗 \itKw')
    % % title('迭代曲线')
    % %
    % % figure(7);
    % % plot(NSGA_II_Result.curve.min(1, :), 'blue', 'LineStyle', '--','LineWidth', 1.3)
    % % hold on
    % % plot(INSGA_II_2_Result.curve.min(1, :), 'black','LineStyle', '-.', 'LineWidth', 1.3)
    % % hold on
    % % grid on; box on;
    % % legend({'NSGA-II', 'I-NSGA-II-QL'}, 'FontName', 'Times New Roman')
    % % xlabel('迭代次数')
    % % ylabel('最大完工时间 \itKw')
    % % title('迭代曲线')
    % %
    % % figure(8);
    % % plot(NSGA_II_Result.curve.min(2, :), 'blue', 'LineStyle', '--','LineWidth', 1.3)
    % % hold on
    % % plot(INSGA_II_2_Result.curve.min(2, :), 'black','LineStyle', '-.', 'LineWidth', 1.3)
    % % hold on
    % % grid on; box on;
    % % legend({'NSGA-II', 'I-NSGA-II-QL'}, 'FontName', 'Times New Roman')
    % % xlabel('迭代次数')
    % % ylabel('总能耗 \itKw')
    % % title('迭代曲线')
    % %
    % % %时间使用反向
    % % figure(9);
    % % plot(NSGA_II_Result.curve.min(1, :), 'blue', 'LineStyle', '--','LineWidth', 1.3)
    % % hold on
    % % plot(INSGA_II_0_Result.curve.min(1, :), 'red', 'LineStyle', '-','LineWidth', 1.3)
    % % hold on
    % % plot(INSGA_II_1_Result.curve.min(1, :), 'magenta','LineStyle', ':', 'LineWidth', 1.3)
    % % hold on
    % % grid on; box on;
    % % legend({'NSGA-II', 'I-NSGA-II-ML', 'I-NSGA-II-OBL'}, 'FontName', 'Times New Roman')
    % % xlabel('迭代次数')
    % % ylabel('最大完工时间 \ith')
    % % title('迭代曲线')
    % %
    % % %时间使用局部
    % % figure(10);
    % % plot(NSGA_II_Result.curve.min(1, :), 'blue', 'LineStyle', '--','LineWidth', 1.3)
    % % hold on
    % % plot(INSGA_II_0_Result.curve.min(1, :), 'red', 'LineStyle', '-','LineWidth', 1.3)
    % % hold on
    % % plot(INSGA_II_2_Result.curve.min(1, :), 'black', 'LineStyle', '-.','LineWidth', 1.3)
    % % grid on; box on;
    % % legend({'NSGA-II', 'I-NSGA-II-ML', 'I-NSGA-II-QL'}, 'FontName', 'Times New Roman')
    % % xlabel('迭代次数')
    % % ylabel('最大完工时间 \ith')
    % % title('迭代曲线')
    % %
    % % %时间使用反向
    % % figure(11);
    % % plot(NSGA_II_Result.curve.min(2, :), 'blue', 'LineStyle', '--','LineWidth', 1.3)
    % % hold on
    % % plot(INSGA_II_0_Result.curve.min(2, :), 'red','LineStyle', '-', 'LineWidth', 1.3)
    % % hold on
    % % plot(INSGA_II_1_Result.curve.min(2, :), 'magenta', 'LineStyle', ':','LineWidth', 1.3)
    % % hold on
    % % grid on; box on;
    % % legend({'NSGA-II', 'I-NSGA-II-ML', 'I-NSGA-II-OBL'}, 'FontName', 'Times New Roman')
    % % xlabel('迭代次数')
    % % ylabel('总能耗 \itKw')
    % % title('迭代曲线')
    % %
    % % %能耗使用局部
    % % figure(12);
    % % plot(NSGA_II_Result.curve.min(2, :), 'blue', 'LineStyle', '--','LineWidth', 1.3)
    % % hold on
    % % plot(INSGA_II_0_Result.curve.min(2, :), 'red','LineStyle', '-', 'LineWidth', 1.3)
    % % hold on
    % % plot(INSGA_II_2_Result.curve.min(2, :), 'black', 'LineStyle', '-.','LineWidth', 1.3)
    % % grid on; box on;
    % % legend({'NSGA-II', 'I-NSGA-II-ML', 'I-NSGA-II-QL'}, 'FontName', 'Times New Roman')
    % % xlabel('迭代次数')
    % % ylabel('总能耗 \itKw')
    % % title('迭代曲线')
    % % %72，
    % % %          % 保存图形句柄
    % % %         figure_handles = [fig1, fig2];
    % % %         % 在每次迭代之间添加一些延迟，以便观察图形
    % % %         pause(2);
    % %
    % % %
    % % %% Pareto图
    % % figure(20)
    % % scatter(NSGA_II_Result.obj_matrix(:,1), NSGA_II_Result.obj_matrix(:,2), 45, 'blue', 's','filled')
    % % hold on;
    % % scatter(INSGA_II_0_Result.obj_matrix(:,1), INSGA_II_0_Result.obj_matrix(:,2), 50, 'red', 'p', 'filled')
    % % hold on;
    % % scatter(INSGA_II_1_Result.obj_matrix(:,1), INSGA_II_1_Result.obj_matrix(:,2), 45, 'magenta', 'd','filled')
    % % hold on;
    % % scatter(INSGA_II_2_Result.obj_matrix(:,1), INSGA_II_2_Result.obj_matrix(:,2), 45, 'black', 'h', 'filled')
    % % grid on; box on;
    % % legend({'NSGA-II', 'I-NSGA-II-ML', 'I-NSGA-II-OBL', 'I-NSGA-II-QL'}, 'FontName', 'Times New Roman',  'Location', 'bestoutside')
    % % xlabel('最大完工时间 Cmax')
    % % ylabel('电量消耗 Kw')
    % % title('Pareto解')
    % % box on; grid on
    % % numresult_I = [numresult_I;INSGA_II_0_Result.obj_matrix];
    % % %numresult_NSGA = [numresult_NSGA;NSGA_II_Result.obj_matrix];
    % %     end
    % %     resultCell_I{j} = numresult_I;
    % %     %resultCell_NSGA{j} = numresult_NSGA;
    % % end
    % % beep;
    %
    %
    % %% 多目标评价指标
    % % 规范化（归一化）
    % total_matrix = [INSGA_II_0_Result.obj_matrix; INSGA_II_3_Result.obj_matrix;];
    % max_obj = max(total_matrix, [], 1);
    % min_obj = min(total_matrix, [], 1);
    % % INSGA2-1 规范化
    % INSGA_II_0_obj_normal = (INSGA_II_0_Result.obj_matrix - repmat(min_obj, size(INSGA_II_0_Result.obj_matrix, 1), 1))./...
    %     (repmat(max_obj, size(INSGA_II_0_Result.obj_matrix, 1), 1) - repmat(min_obj, size(INSGA_II_0_Result.obj_matrix, 1),1));
    % INSGA_II_3_obj_normal = (INSGA_II_3_Result.obj_matrix - repmat(min_obj, size(INSGA_II_3_Result.obj_matrix, 1), 1))./...
    %     (repmat(max_obj, size(INSGA_II_3_Result.obj_matrix, 1), 1) - repmat(min_obj, size(INSGA_II_3_Result.obj_matrix, 1),1));
    %
    % %% HV 指标
    % ref_point = [1.1; 1.1]; % 参考点
    % cd('HV\')
    % HV_ =  [test_lebesgue_measure(INSGA_II_0_obj_normal, ref_point), ...
    %     test_lebesgue_measure(INSGA_II_3_obj_normal, ref_point)];
    % cd('..\')
    % fprintf('HV指标: I-NSGA-II-ML: %.6f \n',HV_(1))
    % fprintf('HV指标: NSGA-II-LS: %.6f \n',HV_(2))
    %
    % %% Spacing
    % cd('Spacing\')
    % Spacing_ =  [Spacing(INSGA_II_0_obj_normal), Spacing(INSGA_II_3_obj_normal)];
    % fprintf('Spacing指标: NSGA-II-ML: %.6f || NSGA-II-LS: %.6f\n',Spacing_(1), Spacing_(2))
    % cd('..\')
    % beep;

    %--------------------------------------------------------------------------%%
    %% 甘特图
    % 需展示的方案的索引
    solution_index = 1;     % 对应的第几个解
    figure(16)
    machine_AGV_gantt_chart(INSGA_II_0_Result.machineTable{solution_index}, INSGA_II_0_Result.AGVTable{solution_index}, ...
        INSGA_II_0_Result.chrom(solution_index, :), jobNum, operaNumVec, AGVSpeed)
    xlabel('Time')
    ylabel('Equipment')
    title(['Makespan：' num2str(INSGA_II_0_Result.obj_matrix(solution_index, 1)) ...
        '\ith \rm Total Energy Consumption：' num2str(INSGA_II_0_Result.obj_matrix(solution_index, 2)) '\itKw·h'])
    fig_filename = fullfile(figureSaveFolder, ['Figure3_' num2str(i) '.fig']);
    savefig(fig_filename);
    fig_filename = fullfile(figureSaveFolder, ['Figure3_' num2str(i) '.png']);
    saveas(gcf, fig_filename, 'png');
end
beep;