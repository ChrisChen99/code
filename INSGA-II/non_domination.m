function f=non_domination(chromosome, variables_number, obj_number)
[N,~]=size(chromosome);
%clear m;
front=1;
F(front).f=[];
individual=[];

%����ÿ�������Ӧ������ֵ
for i=1:N
    individual(i).n=0;    %n  ��Ӧ����i��֧��ĸ�������  ����֧���
    individual(i).p=[];   %p  ��Ӧ����i֧��ĸ��弯��
    for j=1:N
        dom_less=0;
        dom_equal=0;
        dom_more=0;
        for k=1:obj_number
            if chromosome(i,variables_number+k)<chromosome(j,variables_number+k)
                dom_less=dom_less+1;
            elseif chromosome(i,variables_number+k)==chromosome(j,variables_number+k)
                dom_equal=dom_equal+1;
            else
                dom_more=dom_more+1;
            end
        end
        if dom_less==0&&dom_equal~=obj_number  %(�Ż�Ŀ��Ϊ��Сֵ)����i��j������Ŀ��ֵ��С�ڵ���
                                               %����ȫ�����ڣ�i,˵��i��j֧�䡣��i���ڷ�֧��ȼ�1
            individual(i).n=individual(i).n+1;
        elseif dom_more==0&&dom_equal~=obj_number  %(�Ż�Ŀ��Ϊ��Сֵ)����i��j��Ŀ��ֵ�����ڴ���
                                                    %����ȫ�����ڣ�i,˵��i֧��j�����j���뵽i��֧�伯��
            individual(i).p=[individual(i).p j];
        end
    end
    if individual(i).n==0   %����i�ķ�֧�����С ,����֧��ȼ���ߣ����ڵ�ǰ�����е����Ž⣬����Ⱦɫ�м���������Ϣ
        chromosome(i,variables_number+obj_number+1)=1;
        F(front).f=[F(front).f i];  %�ȼ�Ϊ1�ķ�֧��⼯
    end
end
%���ϵĴ����������Ⱦɫ�壬�ҵ��˵ȼ���ߵķ�֧�伯
%                              ÿ������ı�֧����
%                              ÿ�������֧�伯

%%%����Ĵ��뽫����зּ�
while ~isempty(F(front).f)
    Q=[];   %�����һ��front����
    for i=1:length(F(front).f)    %ѭ����ǰ��֧��⼯�еĸ���
        if ~isempty(individual(F(front).f(i)).p)    %����i���Լ���֧��Ľ⼯
            for j=1:length(individual(F(front).f(i)).p)   %ѭ������i��֧��⼯�еĸ���
                individual(individual(F(front).f(i)).p(j)).n=...
                    (individual(individual(F(front).f(i)).p(j)).n)-1;
                if individual(individual(F(front).f(i)).p(j)).n==0
                    chromosome(individual(F(front).f(i)).p(j),...
                        variables_number+obj_number+1)=front+1;
                    Q=[Q individual(F(front).f(i)).p(j)];
                end
            end
        end
    end
    front=front+1;
    F(front).f=Q;
end

[~,index_of_fronts]=sort(chromosome(:,variables_number+obj_number+1));   %d�Ը���Ĵ���������н�����������
for i=1:length(index_of_fronts)
    sorted_based_on_front(i,:)=chromosome(index_of_fronts(i),:);   %sorted_based_on_front�ǰ��յȼ������ľ���
end


%%%Crowding distance ����ÿ�������ӵ����
current_index=0;

for front=1:(length(F)-1)   %%����length(F)-1������ȼ�
    distance=0;
    y=[];
    previous_index=current_index+1;
    for i=1:length(F(front).f)
        y(i,:)=sorted_based_on_front(current_index+i,:);   %y���ŷŵ������򼯺���front�ľ���
    end
    current_index=current_index+i;
    sorted_based_on_objective=[];
    for i=1:obj_number
        [sorted_based_on_objective,index_of_objectives]=sort(y(:,variables_number+i));
        sorted_based_on_objective=[];
        for j=1:length(index_of_objectives)
            sorted_based_on_objective(j,:)=y(index_of_objectives(j),:);
        end
        f_max=sorted_based_on_objective(length(index_of_objectives),variables_number+i);
        f_min=sorted_based_on_objective(1,variables_number+i);
        y(index_of_objectives(length(index_of_objectives)),variables_number+obj_number+1+i)...
            =Inf;
        y(index_of_objectives(1),variables_number+obj_number+1+i)=Inf;
        for j=2:length(index_of_objectives)-1
            next_obj=sorted_based_on_objective(j+1,variables_number+i);
            previous_obj=sorted_based_on_objective(j-1,variables_number+i);
            if (f_max-f_min)==0
                y(index_of_objectives(j),variables_number+obj_number+1+i)=Inf;
            else
                y(index_of_objectives(j),variables_number+obj_number+1+i)=...
                    (next_obj-previous_obj)/(f_max-f_min);
            end
        end
    end
    distance=[];
    distance(:,1)=zeros(length(F(front).f),1);
    for i=1:obj_number
        distance(:,1)=distance(:,1)+y(:,obj_number+variables_number+1+i);
    end
    y(:,obj_number+variables_number+2)=distance;
    y=y(:,1:obj_number+variables_number+2);

    % y����crowed_distanceӵ���ȵĽ�������
    [~,distance_idx]=sort(y(:,obj_number+variables_number+2),'descend');
    y=y(distance_idx,:);

    z(previous_index:current_index,:)=y;
end
f=z();
%%�õ����ǰ����ȼ���ӵ���ȵ���Ⱥ���󣬲����Ѿ����յȼ�����
