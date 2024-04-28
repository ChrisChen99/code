function total_pop = improved_elitism(pop,reverse_percent,jobNum, jobInfo, operaVec, machineNum, AGVNum, AGVSpeed, candidateMachine, ...
    distance_matrix, machineEnergy, AGVEnergy, AGVEG_MAX, AGVEG_MIN, eChargeSpeed, obj_num)
N = size(pop, 1);           % 种群规模
N_reverse = round(N * reverse_percent);   

elitism_pop = pop(1:N_reverse,:);
AGVSpeedNum = length(AGVSpeed); % AGV速度挡位

dim = 5 * sum(operaVec); 
reversed_pop = [];
strong_pop = [];
strongnum = 0;

for i = 1: N_reverse
  
    new_p = reverse_pop(elitism_pop(i, 1: dim), jobNum, operaVec, candidateMachine, AGVNum, AGVSpeedNum);
    old_obj = elitism_pop(i, dim + 1: dim + obj_num);
    new_obj = fitness(new_p, jobNum, jobInfo, operaVec, machineNum, AGVNum, AGVSpeed, candidateMachine, ...
        distance_matrix, machineEnergy, AGVEnergy, AGVEG_MAX, AGVEG_MIN, eChargeSpeed);
    new_obj = new_obj{1};
    if ~ weakly_dominates(old_obj, new_obj)
           reversed_pop = [reversed_pop; [new_p new_obj]];
    end

end
num_individuals_in_reverse_pop = size(reversed_pop, 1);
total_pop = [pop(:, 1: dim + obj_num); reversed_pop];
total_pop = non_domination(total_pop, dim, obj_num);
total_pop = total_pop(1: N, :);
end