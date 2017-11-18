function project_algo1
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
gene = input('number of genes to be generated : ');             %input of number of genes
G = zeros(gene,instances+1);

%randomly make genes and store with objective function value
for ith=1:1:gene
    for j=1:1:instances
        G(ith,j)= random1tok();
    end
    for l= 1:1:con1
        G(ith,C2(l,2)) = G(ith,C2(l,1));
    end
    
    G(ith,instances+1) = objective(G(ith,1:instances));
    %G
end

G = sortrows(G,instances+1);                %sort G matrix wrt objective fn value

%randomly chnge 10% of the value and sort and keep best two everytime
for i=1:1:100
    for j=1:1:gene-2
        ran = random1toinstances();
        G(j,1:instances) = G(j+2,1:instances);
        for l = 1:1:ceil(instances/10)
            G(j,ran(l)) = random1tok();
        end
        for l= 1:1:con1
            G(j,C2(l,2)) = G(j,C2(l,1));
        end
        G(j,instances+1) = objective(G(j,1:instances));
    end
    G = sortrows(G,instances+1);
end

% make 6 genes out of consecutive 2 genes ans select best 2 out of 6
%2 same as selected + 2 with keeping same value in both genes as it is
%+ 2 with keeping different one as it is and change the same one randomly
for i=1:1:50
    for j=1:1:floor(gene/2)
        new = zeros(6,instances+1);
        new(1,1:instances+1) = G(2*(j-1)+1,1:instances+1);
        new(2,1:instances+1) = G(2*(j-1)+2,1:instances+1);
        for l=1:1:instances
            if new(1,l)== new(2,l)
                new(3,l) = new(1,l);
                new(4,l) = new(2,l);
                if rand()<0.5
                    new(5,l) = new(1,l);
                    new(6,l) = new(2,l);
                else
                    new(5,l) = random1tok();
                    new(6,l) = random1tok();
                end
            else
                new(5,l) = new(1,l);
                new(6,l) = new(2,l);
                if rand()<0.5
                    new(3,l) = new(1,l);
                    new(4,l) = new(2,l);
                else
                    new(3,l) = new(2,l);
                    new(4,l) = new(1,l);
                end
            end
        end
        for k=1:1:4
            for l= 1:1:con1
                new(k+2,C2(l,2)) = new(k+2,C2(l,1));
            end
        end
        new(3,instances+1) = objective(new(3,1:instances));
        new(4,instances+1) = objective(new(4,1:instances));
        new(5,instances+1) = objective(new(5,1:instances));
        new(6,instances+1) = objective(new(6,1:instances));
        new = sortrows(new,instances+1);
        G(2*(j-1)+1,1:instances+1) = new(5,1:instances+1);
        G(2*(j-1)+2,1:instances+1) = new(6,1:instances+1);
    end
    G = sortrows(G,instances+1);
end

group = G(gene,:);
group1 = zeros(kclust,1);
element = zeros(kclust,1);
for i = 1:1:instances
    cluster = group(1,i);
    element(cluster,1) = element(cluster,1) + 1;
    group1(cluster,element(cluster,1)) = i;
end
group1

cen = zeros(var,kclust);
    for m = 1:1:instances
        g = group(1,m);
        cen(:,g) = cen(:,g) + A(:,m);
        
    end
    for m = 1:1:kclust
        
            cen(:,m) = cen(:,m)/element(m,1);
       
    end
    
    Dis = 0;
    for m = 1:1:instances
       distance1 = dista(A(:,m),cen(:,group(1,m)));
       Dis = Dis + distance1;
    end
    Dis/100                                             % finding Distance

    function dist = dista(mat1,mat2)
        % function to find distance between two data point
        sum = 0;
        for i = 1:1:var
            sum = sum + ((mat1(i)-mat2(i))*(mat1(i)-mat2(i)));
        end
        dist = sqrt(sum);
    end

    function random = random1tok()      %function to generate randomly 1 to k   
        U = kclust*rand();
        randoms = ceil(U);
        if U == 0
            random = 1;
        else
            random = randoms;
        end
    end

    function random = random1toinstances()   % 10% of instances randomly to change
        random = zeros(1,ceil(instances/10));
        for i= 1:1:ceil(instances/10)
            U = instances*rand();
            randoms = ceil(U);
            if U == 0
                random(i) = 1;
            else
                random(i) = randoms;
            end
        end
    end

    function p = P(mat)                      % function to make P matrix
        p = zeros(instances,kclust);
        for i= 1:1:instances
            p(i,mat(i)) = 1;
        end
    end

    function x=X(mat)                        % function to make X matrix
        B = (P(mat)'*P(mat))^(-1/2);
        x = P(mat)*B;
    end

    function obj = objective(mat)            %function to find objective function value
        cluster = X(mat);
        opt = cluster' * A' * A * cluster;
        obj = trace(opt);
    end
end