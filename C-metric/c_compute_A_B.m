function rate=c_compute_A_B(A,B)  %B�б�A֧��İٷֱ�
NB=size(B,1);
NA=size(A,1);
num=0;
for i=1:NB
    for j=1:NA
        if dominates(A(j,:),B(i,:))
            num=num+1;
            break
        end
    end
end
rate=num/NB;
end