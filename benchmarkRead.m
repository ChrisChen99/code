function benchmarkRead(pth)
[fileID, errmsg] = fopen(pth, 'r');

%% ��ʧ��
if isequal(fileID,-1)
    disp(errmsg);
    return;
end

%% �򿪳ɹ�
firstLine=strip(fgetl(fileID));
firstLineVec=split(firstLine);

%% ��ȡ������������Ŀ
jobNum=str2double(firstLineVec{1});
machineNum=str2double(firstLineVec{2});

%% ¼��ȫ��������Ϣ
operaNumVec=[];                % operaNumVec��¼ÿ��Job�Ĺ�����Ŀ
jobInfo=cell(1,jobNum);        
for i=1:jobNum %
    lineInfo=split(strip(fgetl(fileID)));%��ȡ��һ��
    operaNum=str2double(lineInfo{1});
    operaNumVec=[operaNumVec,operaNum];

    %��ΰѱ�׼����ת��Ϊ�׶�����
    countOs=2;
    for j=1:operaNum
        operaInfo=str2double(lineInfo{countOs});
        processVec=ones(1,machineNum)*Inf;
        for k=1:operaInfo
            countOs=countOs+1;
            machine=str2double(lineInfo{countOs});
            countOs=countOs+1;
            processTime=str2double(lineInfo{countOs});
            processVec(machine)=processTime;
        end
        jobInfo{i}=[jobInfo{i};processVec];
        countOs=countOs+1;
    end
end

%% �ر�
status=fclose(fileID);

%% �ر��ļ�ʧ��
if isequal(status,-1)
    disp('close failed');
    return;
end

%% ��ȡjobInfo����Ϣ��ÿ��Job�ĺ�ѡ����
candidateMachine=[];
for i=1:length(jobInfo)
    for j=1:size(jobInfo{i},1)
        candidateMachine{i,j}=find(jobInfo{i}(j,:)<Inf);
    end
end

%% ����
save('data.mat', 'jobInfo', 'candidateMachine', 'machineNum', 'jobNum', 'operaNumVec');
end