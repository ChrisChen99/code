clc
clear
close all
% ָ�����ڱ���.figͼ�ε��ļ���
figureSaveFolder = 'figures';

if ~exist(figureSaveFolder, 'dir')
    mkdir(figureSaveFolder);
end
numrun = 2;
% param_settings = [
% %     0.7, 0.1, 0.1,0.4,0.7,0.1,0.7;%3
% %     0.7, 0.2, 0.2,0.5,0.8,0.2,0.8;%1(����)
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
%     p_cross = param_settings(j, 1);             % �������
%     p_mutation = param_settings(j, 2);          % �������
%     min_pr = param_settings(j, 3);     % �������ɱ���
%     max_pr = param_settings(j, 4);     % �������ɱ���
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
    disp(['��ǰ�ǵ� ' num2str(i) ' �ζ�������']);
    %
    %         disp(['��' num2str(j) '������ĵ�' num2str(i) '�ζ�������']);

    p_cross = 0.8;              % �������
    p_mutation = 0.2;           % �������
    min_pr =0.1;              % �������ɱ�������
    max_pr = 0.3;              % �������ɱ�������
    epsilon = 0.8;
    alpha = 0.1;
    gamma = 0.8;

    %% �����׼������ʱ������
    pth = 'fjsp\Brandimarte_Data\Mk02.fjs';%2,4,8
    %pth = 'fjsp\Dauzere_Data\18a.fjs';%4,8
    %pth = 'fjsp\Barnes\mt10x.fjs';%4,8
    benchmarkRead(pth);
    load data.mat

    %% �����/װж��վ�������
    %
    % ʹ�����꣨x, y�����ɾ�������
    distance_from_xy(machineNum);
    distance_matrix_excel = xlsread('��������.xlsx', 'װжվ����������');
    % װ��վ��ÿ̨�����ġ�����
    distance_matrix.load_to_machine = distance_matrix_excel(1, :);
    distance_matrix.load_to_machine = distance_matrix.load_to_machine(1: machineNum);
    % ж��վ��ÿ̨�����ľ���
    distance_matrix.machine_to_unload = distance_matrix_excel(2, :);
    distance_matrix.machine_to_unload = distance_matrix.machine_to_unload(1: machineNum);
    % ÿ̨������ľ���
    distance_matrix.machine_to_machine = xlsread('��������.xlsx', '��������������');
    distance_matrix.machine_to_machine = distance_matrix.machine_to_machine(1: machineNum, 1: machineNum);
    % װ ж��վ֮��ľ���
    distance_matrix.load_to_unload = xlsread('��������.xlsx', 'װ��վ��ж��վ����');
    % �ܺ�
    machineEnergy.work = xlsread('��������.xlsx', '�����ӹ��ܺ�');
    machineEnergy.free = xlsread('��������.xlsx', '���������ܺ�');

    %% AGV�������
    AGVNum = 3;
    %AGVSpeed = [0.5,0.75,1];
    AGVSpeed = [0.5,0.75,1];
    % AGV�ܺ�
    AGVEnergy_excel = xlsread('AGV����.xlsx', 'AGV�ܺ�');
    AGVEnergy.free = AGVEnergy_excel(1, :);
    AGVEnergy.load = AGVEnergy_excel(2, :);

    clear distance_matrix_excel AGVEnergy_excel

    %% ���׮����
    % [1]. ���׮������ ж��վ ==> ���׮λ����Ϣ��ж��վ��ͬ
    % AGV����� Kw*h
    AGVEG_MAX = 100;
    % AGV�����ֵ Kw*h �L�L�L�L�L
    % mk01 20kW*h
    % mk02 19.2kW*h
    % mk03 kW*h
    % mk04 24kW*h
    % mk08 24kW*h
    AGVEG_MIN = 19.2;
    % ��֤�����ֵ�Ƿ����
    distance_MAX = max([max(distance_matrix.machine_to_machine) ...
        distance_matrix.load_to_machine ...
        distance_matrix.machine_to_unload...
        distance_matrix.load_to_unload]);
    % ���������ֵ
    % 1.1��ʾ������ 10% ����
    check_MIN = 1.0 * distance_MAX / AGVSpeed(end) * (AGVEnergy.free(end) + AGVEnergy.load(end));
    disp(['������ֵ > ' num2str(check_MIN)])
    if check_MIN > AGVEG_MIN
        error('AGV�����ֵ ==> ���ô���')
    end
    % ����ٶ� Kw
    eChargeSpeed = 20;

    %% NSAG-II �㷨
    % jobNum            ��������
    % jobInfo           �����ӹ�ʱ����Ϣ
    % operaNumVec       ÿ�������Ĺ�������
    % machineNum        ��������
    % AGVNum            AGV����
    % AGVSpeed          AGV�ٶȣ���ͬ��λ��
    % candidateMachine  ÿ�������ĺ�ѡ�ӹ�����
    % distance_matrix   �������
    % machineEnergy     �����ܺ�
    % AGVEnergy         AGV�ܺģ���ͬ��λ��
    % AGVEG_MAX         AGV�����
    % AGVEG_MIN         AGV�����ֵ�����ڴ˵���������г�磩
    % eChargeSpeed      AGV�������
    pop = 100;
    max_gen = 200;
    speedNum = length(AGVSpeed);    % AGV�ٶ�λ��Ŀ
    operaNum = sum(operaNumVec);       % ��������
    %%��Ⱥ��ʼ��
    chrom = init(pop, jobNum, operaNumVec, candidateMachine, AGVNum, speedNum);

    % �㷨�Ա�
    %%%��1��NSGA
    % cd('NSGA-II\')
    % NSGA_II_Result = NSGA2(p_cross,p_mutation,pop,chrom,max_gen,jobNum, jobInfo, operaNumVec, machineNum, AGVNum, AGVSpeed, candidateMachine, distance_matrix, machineEnergy, ...
    %     AGVEnergy, AGVEG_MAX, AGVEG_MIN, eChargeSpeed);
    % cd('..\')
    %
    % % ��3������
    % cd('INSGA-II\')
    % INSGA_II_1_Result = INSGA_II(p_cross,p_mutation,min_pr,max_pr,epsilon,alpha,gamma,pop,chrom,max_gen,jobNum, jobInfo, operaNumVec, machineNum, AGVNum, AGVSpeed, candidateMachine, distance_matrix, machineEnergy, ...
    %     AGVEnergy, AGVEG_MAX, AGVEG_MIN, eChargeSpeed, 'i-elitism');
    % cd('..\')

    % % %��1��i-elitism + ����Q-Learning��VNS
    % cd('INSGA-II\')
    % INSGA_II_3_Result = INSGA_II(p_cross,p_mutation,min_pr,max_pr,epsilon,alpha,gamma,pop,chrom,max_gen,jobNum, jobInfo, operaNumVec, machineNum, AGVNum, AGVSpeed, candidateMachine, distance_matrix, machineEnergy, ...
    %     AGVEnergy, AGVEG_MAX, AGVEG_MIN, eChargeSpeed, 'NOQ_VNS+i-elitism');
    % cd('..\')
    %
    %��1��i-elitism + VNS
    cd('INSGA-II\')
    INSGA_II_0_Result = INSGA_II(p_cross,p_mutation,min_pr,max_pr,epsilon,alpha,gamma,pop,chrom,max_gen,jobNum, jobInfo, operaNumVec, machineNum, AGVNum, AGVSpeed, candidateMachine, distance_matrix, machineEnergy, ...
        AGVEnergy, AGVEG_MAX, AGVEG_MIN, eChargeSpeed, 'VNS+i-elitism');
    cd('..\')
    %
    % %%��2��VNS
    % cd('INSGA-II\')
    % INSGA_II_2_Result = INSGA_II(p_cross,p_mutation,min_pr,max_pr,epsilon,alpha,gamma,pop,chrom,max_gen,jobNum, jobInfo, operaNumVec, machineNum, AGVNum, AGVSpeed, candidateMachine, distance_matrix, machineEnergy, ...
    %     AGVEnergy, AGVEG_MAX, AGVEG_MIN, eChargeSpeed, 'VNS');
    % cd('..\')


    %
    % figure(2);
    % plot(NSGA_II_Result.curve.min(1, :), 'Color', [0, 0.447, 0.741], 'LineStyle', '-', 'LineWidth', 1.5) % ��ɫ
    % hold on
    % plot(INSGA_II_2_Result.curve.min(1, :), 'Color', [0.85, 0.325, 0.098], 'LineStyle', '-', 'LineWidth', 1.5) % ��ɫ
    % hold on
    % grid on; box on;
    % legend({'NSGA-II', 'I-NSGA-II-ML'}, 'FontName', 'Times New Roman')
    % xlabel('��������', 'FontName', 'Times New Roman')
    % ylabel('����깤ʱ�� \ith', 'FontName', 'Times New Roman')
    % set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)
    % set(gcf, 'Position', [100, 100, 800, 600])  % ����ͼ�Ĵ�С


    % % ����һ��ͼ�δ���
    % figure(1);
    % % ��������
    % x1 = NSGA_II_Result.curve.min(1, :);
    % x2 = INSGA_II_2_Result.curve.min(1, :);
    % % ������ɫ
    % colors = {[0, 0.447, 0.741], [0.85, 0.325, 0.098]}; % ��ɫ�ͺ�ɫ
    % % ������������
    % plot(x1, 'Color', colors{1}, 'LineStyle', '-', 'LineWidth', 1.5);
    % hold on;
    % plot(x2, 'Color', colors{2}, 'LineStyle', '-', 'LineWidth', 1.5);
    % % ���ÿ10�����ӵ�5����ʼ
    % indices = 5:10:numel(x1);
    % for j = 1:length(indices)
    %     plot(indices(j), x1(indices(j)), 's', 'MarkerSize', 10, 'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', 'k');
    %     plot(indices(j), x2(indices(j)), 's', 'MarkerSize', 10, 'MarkerFaceColor', colors{2}, 'MarkerEdgeColor', 'k');
    % end
    % grid on;
    % box on;
    % legend({'','','NSGA-II', 'I-NSGA-II-ML'}, 'FontName', 'Times New Roman')
    % % ����ͼ��
    % xlabel('��������', 'FontName', 'Times New Roman');
    % ylabel('����깤ʱ�� \ith', 'FontName', 'Times New Roman');
    % set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
    % set(gcf, 'Position', [100, 100, 800, 600]);
    % set(gca,'Box','off');
    %
    % % %% ��������
    % %
    % % %         fig1 = figure('Name', ['ʱ��-��' num2str(j) '�Ĳ�����' num2str(i) '���������']);
    % % %         % ��������ӻ�ͼ�����ݲ���
    % % figure(1);
    % % plot(NSGA_II_Result.curve.min(1, :), 'blue', 'LineStyle', '--','LineWidth', 1.3)
    % % hold on
    % % plot(INSGA_II_0_Result.curve.min(1, :), 'red', 'LineStyle', '-','LineWidth', 1.3)
    % % hold on
    % % %
    % % %         %����
    % % %         legend({'NSGA-II', 'I-NSGA-II-ML'}, 'FontName', 'Times New Roman')
    % %
    % % plot(INSGA_II_1_Result.curve.min(1, :), 'magenta','LineStyle', ':', 'LineWidth', 1.3)
    % % hold on
    % % plot(INSGA_II_2_Result.curve.min(1, :), 'black', 'LineStyle', '-.','LineWidth', 1.3)
    % % grid on; box on;
    % % legend({'NSGA-II', 'I-NSGA-II-ML', 'I-NSGA-II-OBL', 'I-NSGA-II-QL'}, 'FontName', 'Times New Roman')
    % % xlabel('��������')
    % % ylabel('����깤ʱ�� \ith')
    % % title('��������')
    % %
    % % %
    % % %fig2 = figure('Name', ['�ܺ�-��' num2str(j) '�����ĵ�' num2str(i) '���������']);
    % % % ��������ӻ�ͼ�����ݲ���
    % % figure(2);
    % % plot(NSGA_II_Result.curve.min(2, :), 'blue', 'LineStyle', '--','LineWidth', 1.3)
    % % hold on
    % % plot(INSGA_II_0_Result.curve.min(2, :), 'red','LineStyle', '-', 'LineWidth', 1.3)
    % % hold on
    % % %%����
    % % %legend({'NSGA-II', 'I-NSGA-II-ML'}, 'FontName', 'Times New Roman')
    % %
    % % plot(INSGA_II_1_Result.curve.min(2, :), 'magenta', 'LineStyle', ':','LineWidth', 1.3)
    % % hold on
    % % plot(INSGA_II_2_Result.curve.min(2, :), 'black', 'LineStyle', '-.','LineWidth', 1.3)
    % % grid on; box on;
    % % legend({'NSGA-II', 'I-NSGA-II-ML', 'I-NSGA-II-OBL', 'I-NSGA-II-QL'}, 'FontName', 'Times New Roman')
    % % xlabel('��������')
    % % ylabel('���ܺ� \itKw')
    % % title('��������')
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
    % % xlabel('��������')
    % % ylabel('����깤ʱ�� \itKw')
    % % title('��������')
    % %
    % % figure(4);
    % % plot(NSGA_II_Result.curve.min(2, :), 'blue', 'LineStyle', '--','LineWidth', 1.3)
    % % hold on
    % % plot(INSGA_II_0_Result.curve.min(2, :), 'red','LineStyle', '-', 'LineWidth', 1.3)
    % % hold on
    % % grid on; box on;
    % % legend({'NSGA-II', 'I-NSGA-II-ML'}, 'FontName', 'Times New Roman')
    % % xlabel('��������')
    % % ylabel('���ܺ� \itKw')
    % % title('��������')
    % %
    % % figure(5);
    % % plot(NSGA_II_Result.curve.min(1, :), 'blue', 'LineStyle', '--','LineWidth', 1.3)
    % % hold on
    % % plot(INSGA_II_1_Result.curve.min(1, :), 'magenta','LineStyle', ':', 'LineWidth', 1.3)
    % % hold on
    % % grid on; box on;
    % % legend({'NSGA-II', 'I-NSGA-II-OBL'}, 'FontName', 'Times New Roman')
    % % xlabel('��������')
    % % ylabel('����깤ʱ�� \itKw')
    % % title('��������')
    %
    % % figure(5);
    % % plot(NSGA_II_Result.curve.min(1, :), 'b--', 'LineWidth', 1.3)
    % % hold on
    % % plot(INSGA_II_1_Result.curve.min(1, :), 'm:', 'LineWidth', 1.3)
    % % grid on; box on;
    % % legend({'NSGA-II', 'I-NSGA-II-OBL'}, 'FontName', 'Times New Roman', 'FontSize', 12)
    % % xlabel('��������', 'FontName',  'FontSize', 12)
    % % ylabel('����깤ʱ�� \itKw', 'FontName', 'FontSize', 12)
    % % title('��������', 'FontName', 'FontSize', 14)
    % % set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)
    % % set(gcf, 'Position', [100, 100, 800, 600])  % ����ͼ�Ĵ�С
    % %
    % % figure(5);
    % % % ʹ�ò�ͬ����ɫӳ��
    % % cmap = colormap(parula(2)); % ����ʹ����Parula��ɫ�壬����2��ʾ������ɫ
    % % colormap(cmap);
    % % plot(NSGA_II_Result.curve.min(1, :), 'Color', cmap(2, :), 'LineWidth', 1.3) % ָ����ɫΪcmap(1, :)
    % % hold on
    % % plot(INSGA_II_1_Result.curve.min(1, :), 'Color', cmap(1, :), 'LineWidth', 1.3) % ָ����ɫΪcmap(2, :)
    % % % ������ݵ���
    % % scatter(1:length(NSGA_II_Result.curve.min(1, :)), NSGA_II_Result.curve.min(1, :), 40, cmap(2, :), 'v', 'filled')
    % % scatter(1:length(INSGA_II_1_Result.curve.min(1, :)), INSGA_II_1_Result.curve.min(1, :), 40, cmap(1, :), 'v', 'filled')
    % % grid on; box on;
    % % legend({'NSGA-II', 'I-NSGA-II-OBL'}, 'FontName', 'Times New Roman', 'FontSize', 12)
    % % xlabel('��������', 'FontName', '����', 'FontSize', 12)
    % % ylabel('����깤ʱ�� \ith', 'FontName', '����', 'FontSize', 12)
    % % title('��������', 'FontName', '����', 'FontSize', 14)
    % % set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)
    % % set(gcf, 'Position', [100, 100, 800, 600])  % ����ͼ�Ĵ�С
    % %
    % % figure(6);
    % % % ʹ�ò�ͬ����ɫӳ��
    % % cmap = colormap(parula(2)); % ����ʹ����Parula��ɫ�壬����2��ʾ������ɫ
    % % colormap(cmap);
    % % plot(NSGA_II_Result.curve.min(2, :), 'Color', cmap(2, :), 'LineWidth', 1.3) % ָ����ɫΪcmap(1, :)
    % % hold on
    % % plot(INSGA_II_1_Result.curve.min(2, :), 'Color', cmap(1, :), 'LineWidth', 1.3) % ָ����ɫΪcmap(2, :)
    % % % ������ݵ���
    % % scatter(1:length(NSGA_II_Result.curve.min(2, :)), NSGA_II_Result.curve.min(2, :), 40, cmap(2, :), 'v', 'filled')
    % % scatter(1:length(INSGA_II_1_Result.curve.min(2, :)), INSGA_II_1_Result.curve.min(2, :), 40, cmap(1, :), 'v', 'filled')
    % % grid on; box on;
    % % legend({'NSGA-II', 'I-NSGA-II-OBL'}, 'FontName', 'Times New Roman', 'FontSize', 12)
    % % xlabel('��������', 'FontName', '����', 'FontSize', 12)
    % % ylabel('���ܺ� \itKw', 'FontName', '����', 'FontSize', 12)
    % % title('��������', 'FontName', '����', 'FontSize', 14)
    % % set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)
    % % set(gcf, 'Position', [100, 100, 800, 600])  % ����ͼ�Ĵ�С
    % % figure(6);
    % % plot(NSGA_II_Result.curve.min(2, :), 'blue', 'LineStyle', '--','LineWidth', 1.3)
    % % hold on
    % % plot(INSGA_II_1_Result.curve.min(2, :), 'magenta','LineStyle', ':', 'LineWidth', 1.3)
    % % hold on
    % % grid on; box on;
    % % legend({'NSGA-II', 'I-NSGA-II-OBL'}, 'FontName', 'Times New Roman')
    % % xlabel('��������')
    % % ylabel('���ܺ� \itKw')
    % % title('��������')
    % %
    % % figure(7);
    % % plot(NSGA_II_Result.curve.min(1, :), 'blue', 'LineStyle', '--','LineWidth', 1.3)
    % % hold on
    % % plot(INSGA_II_2_Result.curve.min(1, :), 'black','LineStyle', '-.', 'LineWidth', 1.3)
    % % hold on
    % % grid on; box on;
    % % legend({'NSGA-II', 'I-NSGA-II-QL'}, 'FontName', 'Times New Roman')
    % % xlabel('��������')
    % % ylabel('����깤ʱ�� \itKw')
    % % title('��������')
    % %
    % % figure(8);
    % % plot(NSGA_II_Result.curve.min(2, :), 'blue', 'LineStyle', '--','LineWidth', 1.3)
    % % hold on
    % % plot(INSGA_II_2_Result.curve.min(2, :), 'black','LineStyle', '-.', 'LineWidth', 1.3)
    % % hold on
    % % grid on; box on;
    % % legend({'NSGA-II', 'I-NSGA-II-QL'}, 'FontName', 'Times New Roman')
    % % xlabel('��������')
    % % ylabel('���ܺ� \itKw')
    % % title('��������')
    % %
    % % %ʱ��ʹ�÷���
    % % figure(9);
    % % plot(NSGA_II_Result.curve.min(1, :), 'blue', 'LineStyle', '--','LineWidth', 1.3)
    % % hold on
    % % plot(INSGA_II_0_Result.curve.min(1, :), 'red', 'LineStyle', '-','LineWidth', 1.3)
    % % hold on
    % % plot(INSGA_II_1_Result.curve.min(1, :), 'magenta','LineStyle', ':', 'LineWidth', 1.3)
    % % hold on
    % % grid on; box on;
    % % legend({'NSGA-II', 'I-NSGA-II-ML', 'I-NSGA-II-OBL'}, 'FontName', 'Times New Roman')
    % % xlabel('��������')
    % % ylabel('����깤ʱ�� \ith')
    % % title('��������')
    % %
    % % %ʱ��ʹ�þֲ�
    % % figure(10);
    % % plot(NSGA_II_Result.curve.min(1, :), 'blue', 'LineStyle', '--','LineWidth', 1.3)
    % % hold on
    % % plot(INSGA_II_0_Result.curve.min(1, :), 'red', 'LineStyle', '-','LineWidth', 1.3)
    % % hold on
    % % plot(INSGA_II_2_Result.curve.min(1, :), 'black', 'LineStyle', '-.','LineWidth', 1.3)
    % % grid on; box on;
    % % legend({'NSGA-II', 'I-NSGA-II-ML', 'I-NSGA-II-QL'}, 'FontName', 'Times New Roman')
    % % xlabel('��������')
    % % ylabel('����깤ʱ�� \ith')
    % % title('��������')
    % %
    % % %ʱ��ʹ�÷���
    % % figure(11);
    % % plot(NSGA_II_Result.curve.min(2, :), 'blue', 'LineStyle', '--','LineWidth', 1.3)
    % % hold on
    % % plot(INSGA_II_0_Result.curve.min(2, :), 'red','LineStyle', '-', 'LineWidth', 1.3)
    % % hold on
    % % plot(INSGA_II_1_Result.curve.min(2, :), 'magenta', 'LineStyle', ':','LineWidth', 1.3)
    % % hold on
    % % grid on; box on;
    % % legend({'NSGA-II', 'I-NSGA-II-ML', 'I-NSGA-II-OBL'}, 'FontName', 'Times New Roman')
    % % xlabel('��������')
    % % ylabel('���ܺ� \itKw')
    % % title('��������')
    % %
    % % %�ܺ�ʹ�þֲ�
    % % figure(12);
    % % plot(NSGA_II_Result.curve.min(2, :), 'blue', 'LineStyle', '--','LineWidth', 1.3)
    % % hold on
    % % plot(INSGA_II_0_Result.curve.min(2, :), 'red','LineStyle', '-', 'LineWidth', 1.3)
    % % hold on
    % % plot(INSGA_II_2_Result.curve.min(2, :), 'black', 'LineStyle', '-.','LineWidth', 1.3)
    % % grid on; box on;
    % % legend({'NSGA-II', 'I-NSGA-II-ML', 'I-NSGA-II-QL'}, 'FontName', 'Times New Roman')
    % % xlabel('��������')
    % % ylabel('���ܺ� \itKw')
    % % title('��������')
    % % %72��
    % % %          % ����ͼ�ξ��
    % % %         figure_handles = [fig1, fig2];
    % % %         % ��ÿ�ε���֮�����һЩ�ӳ٣��Ա�۲�ͼ��
    % % %         pause(2);
    % %
    % % %
    % % %% Paretoͼ
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
    % % xlabel('����깤ʱ�� Cmax')
    % % ylabel('�������� Kw')
    % % title('Pareto��')
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
    % %% ��Ŀ������ָ��
    % % �淶������һ����
    % total_matrix = [INSGA_II_0_Result.obj_matrix; INSGA_II_3_Result.obj_matrix;];
    % max_obj = max(total_matrix, [], 1);
    % min_obj = min(total_matrix, [], 1);
    % % INSGA2-1 �淶��
    % INSGA_II_0_obj_normal = (INSGA_II_0_Result.obj_matrix - repmat(min_obj, size(INSGA_II_0_Result.obj_matrix, 1), 1))./...
    %     (repmat(max_obj, size(INSGA_II_0_Result.obj_matrix, 1), 1) - repmat(min_obj, size(INSGA_II_0_Result.obj_matrix, 1),1));
    % INSGA_II_3_obj_normal = (INSGA_II_3_Result.obj_matrix - repmat(min_obj, size(INSGA_II_3_Result.obj_matrix, 1), 1))./...
    %     (repmat(max_obj, size(INSGA_II_3_Result.obj_matrix, 1), 1) - repmat(min_obj, size(INSGA_II_3_Result.obj_matrix, 1),1));
    %
    % %% HV ָ��
    % ref_point = [1.1; 1.1]; % �ο���
    % cd('HV\')
    % HV_ =  [test_lebesgue_measure(INSGA_II_0_obj_normal, ref_point), ...
    %     test_lebesgue_measure(INSGA_II_3_obj_normal, ref_point)];
    % cd('..\')
    % fprintf('HVָ��: I-NSGA-II-ML: %.6f \n',HV_(1))
    % fprintf('HVָ��: NSGA-II-LS: %.6f \n',HV_(2))
    %
    % %% Spacing
    % cd('Spacing\')
    % Spacing_ =  [Spacing(INSGA_II_0_obj_normal), Spacing(INSGA_II_3_obj_normal)];
    % fprintf('Spacingָ��: NSGA-II-ML: %.6f || NSGA-II-LS: %.6f\n',Spacing_(1), Spacing_(2))
    % cd('..\')
    % beep;

    %--------------------------------------------------------------------------%%
    %% ����ͼ
    % ��չʾ�ķ���������
    solution_index = 1;     % ��Ӧ�ĵڼ�����
    figure(16)
    machine_AGV_gantt_chart(INSGA_II_0_Result.machineTable{solution_index}, INSGA_II_0_Result.AGVTable{solution_index}, ...
        INSGA_II_0_Result.chrom(solution_index, :), jobNum, operaNumVec, AGVSpeed)
    xlabel('Time')
    ylabel('Equipment')
    title(['Makespan��' num2str(INSGA_II_0_Result.obj_matrix(solution_index, 1)) ...
        '\ith \rm Total Energy Consumption��' num2str(INSGA_II_0_Result.obj_matrix(solution_index, 2)) '\itKw��h'])
    fig_filename = fullfile(figureSaveFolder, ['Figure3_' num2str(i) '.fig']);
    savefig(fig_filename);
    fig_filename = fullfile(figureSaveFolder, ['Figure3_' num2str(i) '.png']);
    saveas(gcf, fig_filename, 'png');
end
beep;