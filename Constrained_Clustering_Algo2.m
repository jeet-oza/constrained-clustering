function project_algo2
%Algorithm for constrained clustering
Afile = input('enter the filename for A matrix : ');
Asheet = input('enter the sheet number of A matrix : ');
A = xlsread(Afile,Asheet);                            %input of data matrix
A = A';
[var,instances] = size(A);
b = randperm(instances);

kclust = input('number of cluster : ');            %input number of cluster

Cfile = input('enter the filename for C1 matrix : ');
Csheet = input('enter the sheet number of C1 matrix : ');
C1 = xlsread(Cfile,Csheet);                    %input of constrained matrix
con =1;

%finding no. of constraint and store
for i= 1:1:instances
    for j= 1:1:i
        if C1(i,j) == 1 
            C2(con,1) = i;
            C2(con,2) = j;
            con = con+1;
        end
    end
end

[con1,con2] = size(C2);



%make k centroid for each parameter varying from mean-2*std to mean+2*std
cen = zeros(var,kclust);
k = 4/(kclust-1);
%for n = 1:1:var
    for m =1:1:kclust
        cen(:,m) = A(:,b(m));
    end
%end

%algorithm
% step 1 :group element to nearest centroid
% step 2 : find number of element
% step 3 : find new centroid for each group with new element
% step 4: check terminating condition if terminated leave otherwise go to
% step 1
group = zeros(1,instances);
group1 = zeros(kclust,1);
element = zeros(kclust,1);
Distance = zeros(1,2);
while(1)
    number = zeros(1,kclust);
    Dis = 0;
    for m = 1:1:instances
        distance1 = zeros(1,kclust);
        for l = 1:1:kclust
            distance1(1,l) = dista(A(:,m),cen(:,l));
            
        end
        [C,I] = min(distance1);
        Dis = Dis + C;
        group(1,m) = I;
    end
    Distance(1,1) = Distance(1,2);
    Distance(1,2) = Dis;
    gap = Distance(1,1) - Distance(1,2);
    for m= 1:1:con1
        group(1,C2(m,2)) = group(1,C2(m,1));
    end
    cen = zeros(var,kclust);
    for m = 1:1:instances
        g = group(1,m);
        cen(:,g) = cen(:,g) + A(:,m);
        number(1,g) = number(1,g)+1;
    end
    for m = 1:1:kclust
        if number(1,m) == 0
            cen(:,m) = A(:,b(mod(2*m,instances)));
        else
            cen(:,m) = cen(:,m)/number(1,m);
        end
    end
    if gap == 0
        break
    end
end
%group                   %show which element is in which group

for i = 1:1:instances
    cluster = group(1,i);
    element(cluster,1) = element(cluster,1) + 1;
    group1(cluster,element(cluster,1)) = i;
end

group1(1,:)
group1(2,:)
group1(3,:)
group1(4,:)
group1(5,:)
group1(6,:)
Distance/100
    function dist = dista(mat1,mat2)
        % function to find distance between two data point
        sum = 0;
        for i = 1:1:var
            sum = sum + ((mat1(i)-mat2(i))*(mat1(i)-mat2(i)));
        end
        dist = sqrt(sum);
    end

end