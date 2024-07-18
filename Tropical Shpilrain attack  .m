[key,U,V,A,B,W] = MaxplusGenerateKeyStickels(5, -1000, 1000, 50, -1000,1000);
n = size(A, 1);

X_ind = reshape(1:n^2, n, []).';
Y_ind = reshape(n^2+1:2*n^2, n, []).';
T_ind = reshape(2*n^2+1:3*n^2, n, []).';
R_ind = reshape(3*n^2+1:4*n^2, n, []).';
w1_ind = reshape(4*n^2+1:4*n^2+n^3, n, []).';
w2_ind = reshape(4*n^2+n^3+1:4*n^2+2*n^3, n, []).';
w3_ind = reshape(4*n^2+2*n^3+1:4*n^2+3*n^3, n, []).';
w4_ind = reshape(4*n^2+3*n^3+1:4*n^2+4*n^3, n, []).';
w5_ind = reshape(4*n^2+4*n^3+1:4*n^2+4*n^3+n^4, n*n, []).';

M = 1000000;
Ain = zeros(8*n^3+2*n^4, 4*n^2+4*n^3+n^4);
bin = zeros(8*n^3+2*n^4, 1);

Aeq = zeros(5*n^2, 4*n^2+4*n^3+n^4);
beq = zeros(5*n^2, 1);

Ain_row_counter = 1;
Aeq_row_counter = 1;

%case1: XA=T
for i = 1:n
    for j = 1:n
        %<=
        for k = 1:n
            Ain(Ain_row_counter, X_ind(i, k)) = 1;
            Ain(Ain_row_counter, T_ind(i, j)) = -1;
            bin(Ain_row_counter) = -A(k, j);
            Ain_row_counter = Ain_row_counter + 1;
        end

        %>=
        for k = 1:n
            Ain(Ain_row_counter, X_ind(i, k)) = -1;
            Ain(Ain_row_counter, T_ind(i, j)) = 1;
            Ain(Ain_row_counter, w1_ind(n*(i-1)+j, k)) = M;
            bin(Ain_row_counter) = A(k, j) + M;
            Ain_row_counter = Ain_row_counter + 1;
        end
        Aeq(Aeq_row_counter, w1_ind(n*(i-1)+j, 1:n)) = ones(1, n);
        beq(Aeq_row_counter) = 1;
        Aeq_row_counter = Aeq_row_counter + 1;
    end
end

%Case2: AX=T
for i = 1:n
    for j = 1:n
        %<=
        for k = 1:n
            Ain(Ain_row_counter, X_ind(k, j)) = 1;
            Ain(Ain_row_counter, T_ind(i, j)) = -1;
            bin(Ain_row_counter) = -A(i, k);
            Ain_row_counter = Ain_row_counter + 1;
        end

        %>=
        for k = 1:n
            Ain(Ain_row_counter, X_ind(k, j)) = -1;
            Ain(Ain_row_counter, T_ind(i, j)) = 1;
            Ain(Ain_row_counter, w2_ind(n*(i-1)+j, k)) = M;
            bin(Ain_row_counter) = A(i, k) + M;
            Ain_row_counter = Ain_row_counter + 1;
        end
        Aeq(Aeq_row_counter, w2_ind(n*(i-1)+j, 1:n)) = ones(1, n);
        beq(Aeq_row_counter) = 1;
        Aeq_row_counter = Aeq_row_counter + 1;
    end
end

%case3: YB=R
for i = 1:n
    for j = 1:n
        %<=
        for k = 1:n
            Ain(Ain_row_counter, Y_ind(i, k)) = 1;
            Ain(Ain_row_counter, R_ind(i, j)) = -1;
            bin(Ain_row_counter) = -B(k, j);
            Ain_row_counter = Ain_row_counter + 1;
        end

        %>=
        for k = 1:n
            Ain(Ain_row_counter, Y_ind(i, k)) = -1;
            Ain(Ain_row_counter, R_ind(i, j)) = 1;
            Ain(Ain_row_counter, w3_ind(n*(i-1)+j, k)) = M;
            bin(Ain_row_counter) = B(k, j) + M;
            Ain_row_counter = Ain_row_counter + 1;
        end
        Aeq(Aeq_row_counter, w3_ind(n*(i-1)+j, 1:n)) = ones(1, n);
        beq(Aeq_row_counter) = 1;
        Aeq_row_counter = Aeq_row_counter + 1;
    end
end

%Case4: BY=R
for i = 1:n
    for j = 1:n
        %<=
        for k = 1:n
            Ain(Ain_row_counter, Y_ind(k, j)) = 1;
            Ain(Ain_row_counter, R_ind(i, j)) = -1;
            bin(Ain_row_counter) = -B(i, k);
            Ain_row_counter = Ain_row_counter + 1;
        end

        %>=
        for k = 1:n
            Ain(Ain_row_counter, Y_ind(k, j)) = -1;
            Ain(Ain_row_counter, R_ind(i, j)) = 1;
            Ain(Ain_row_counter, w4_ind(n*(i-1)+j, k)) = M;
            bin(Ain_row_counter) = B(i, k) + M;
            Ain_row_counter = Ain_row_counter + 1;
        end
        Aeq(Aeq_row_counter, w4_ind(n*(i-1)+j, 1:n)) = ones(1, n);
        beq(Aeq_row_counter) = 1;
        Aeq_row_counter = Aeq_row_counter + 1;
    end
end

%Case5:XWY=U
for i = 1:n
    for j = 1:n
        %<=
        for k = 1:n
            for l=1:n
            Ain(Ain_row_counter, X_ind(i, k)) = 1;
            Ain(Ain_row_counter, Y_ind(l, j)) = 1;
            bin(Ain_row_counter) = U(i, j)-W(k,l);
            Ain_row_counter = Ain_row_counter + 1;
            end
        end

        %>=
        for k = 1:n
            for l=1:n
            Ain(Ain_row_counter, X_ind(i, k)) = -1;
            Ain(Ain_row_counter, Y_ind(l, j)) = -1;
            Ain(Ain_row_counter, w5_ind(n*(i-1)+j,n*(k-1)+l)) = M;
            bin(Ain_row_counter) = -U(i, j) + M+W(k,l);
            Ain_row_counter = Ain_row_counter + 1;
            end
        end
        Aeq(Aeq_row_counter, w5_ind(n*(i-1)+j, 1:n^2)) = ones(1, n^2);
        beq(Aeq_row_counter) = 1;
        Aeq_row_counter = Aeq_row_counter + 1;
    end
end

f = zeros(1, 4*n^2+4*n^3+n^4);
intcon = 1:4*n^2+4*n^3+n^4;
lb = -inf(4*n^2+4*n^3+n^4, 1);
lb(4*n^2+1:end) = 0;
ub = inf(4*n^2+4*n^3+n^4, 1);
ub(4*n^2+1:end) = 1;

model.A = sparse([Ain; Aeq]);
model.obj = f;
model.modelsense = 'min';
model.rhs = [bin; beq];
model.lb = lb;
model.ub = ub;
model.vtype = [repmat('C', 1, 4*n^2), repmat('B', 1, 4*n^3+n^4)];
model.sense = [repmat('<', size(Ain, 1), 1); repmat('=', size(Aeq, 1), 1)];

params.outputflag = 1;
params.BranchDir=  0;
params.MIPFocus=  1;
params.TimeLimit =60*1000;
result = gurobi(model, params);

z = result.x;

x = z(1:n^2);
y = z(n^2+1:2*n^2);
temp = MaxplusMulti(reshape(x, n, n)', V);
K_attack = MaxplusMulti(temp, reshape(y, n, n)');
key == K_attack
key
K_attack

