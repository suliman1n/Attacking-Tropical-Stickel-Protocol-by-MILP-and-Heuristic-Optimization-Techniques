n_values = [3, 5,7,8];
D_values = [3,5,10,15,20,25,30,35,40,45,50];
successful_trials_needed = 3; 


average_times = zeros(length(n_values), length(D_values));


for n_index = 1:length(n_values)
    n = n_values(n_index);
    
    for D_index = 1:length(D_values)
        D = D_values(D_index);
        
        
        trial_times = zeros(1, successful_trials_needed);
        successful_trials_count = 0;
        
        
        while successful_trials_count < successful_trials_needed
            
           [key,U,V,A,B,W] = GenerateKeyStickelsmaxmin(n, -1000, 1000, D, -1000,1000);

n = size(A, 1);

X_ind = reshape(1:n^2, n, []).';
Y_ind = reshape(n^2+1:2*n^2, n, []).';
T_ind = reshape(2*n^2+1:3*n^2, n, []).';
R_ind = reshape(3*n^2+1:4*n^2, n, []).';
w1_ind = reshape(4*n^2+1:4*n^2+n^3, n, []).';
w2_ind = reshape(4*n^2+n^3+1:4*n^2+2*n^3, n, []).';
w3_ind = reshape(4*n^2+2*n^3+1:4*n^2+3*n^3, n, []).';
w4_ind = reshape(4*n^2+3*n^3+1:4*n^2+4*n^3, n, []).';
w5_ind = reshape(4*n^2+4*n^3+1:4*n^2+4*n^3+n^4, n^2, []).';
wh_ind=reshape(4*n^2+4*n^3+n^4+1:4*n^2+12*n^3+n^4, 2, []).';
wh_ind5=reshape(4*n^2+12*n^3+n^4+1:4*n^2+12*n^3+4*n^4, 3, []).';

M = 10000;
Ain = zeros(16*n^3+6*n^4, 4*n^2+12*n^3+4*n^4);
bin = zeros(16*n^3+6*n^4, 1);

Aeq = zeros(4*n^3+n^4+5*n^2, 4*n^2+12*n^3+4*n^4);
beq = zeros(4*n^3+n^4+5*n^2, 1);

Ain_row_counter = 1;
Aeq_row_counter = 1;
wh_ind_row_counter=1;
wh_ind5_row_counter=1;
%case1: XA=T
for i = 1:n
    for j = 1:n
        %<=
        for k = 1:n
            %eq1
            Ain(Ain_row_counter, X_ind(i, k)) = 1;
            Ain(Ain_row_counter, T_ind(i, j)) = -1;
            Ain(Ain_row_counter, wh_ind(wh_ind_row_counter,1))=M;
            bin(Ain_row_counter) = M;
            Ain_row_counter = Ain_row_counter + 1;
            %eq2
            Ain(Ain_row_counter, T_ind(i, j)) = -1;
            Ain(Ain_row_counter, wh_ind(wh_ind_row_counter,2))=M;
            bin(Ain_row_counter) = M-A(k, j);
            Ain_row_counter = Ain_row_counter + 1;
            %eq3
            Aeq(Aeq_row_counter,wh_ind(wh_ind_row_counter,1))=1;
            Aeq(Aeq_row_counter,wh_ind(wh_ind_row_counter,2))=1;
            beq(Aeq_row_counter) = 1;
            Aeq_row_counter = Aeq_row_counter + 1;
            wh_ind_row_counter=wh_ind_row_counter+1;
        end

        %>=
        for k = 1:n
            %eq1
            Ain(Ain_row_counter, X_ind(i, k)) = -1;
            Ain(Ain_row_counter, T_ind(i, j)) = 1;
            Ain(Ain_row_counter, w1_ind(n*(i-1)+j, k)) = M;
            bin(Ain_row_counter) = M;
            Ain_row_counter = Ain_row_counter + 1;
            %eq2
            Ain(Ain_row_counter, T_ind(i, j)) = 1;
            Ain(Ain_row_counter, w1_ind(n*(i-1)+j, k)) = M;
            bin(Ain_row_counter) = A(k, j)+M;
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
            %eq1
            Ain(Ain_row_counter, T_ind(i, j)) = -1;
            Ain(Ain_row_counter, wh_ind(wh_ind_row_counter,1))=M;
            bin(Ain_row_counter) = -A(i, k)+M;
            Ain_row_counter = Ain_row_counter + 1;
            %eq2
            Ain(Ain_row_counter, X_ind(k, j)) = 1;
            Ain(Ain_row_counter, T_ind(i, j)) = -1;
            Ain(Ain_row_counter, wh_ind(wh_ind_row_counter,2))=M;
            bin(Ain_row_counter) = M;
            Ain_row_counter = Ain_row_counter + 1;
            %eq3
            Aeq(Aeq_row_counter,wh_ind(wh_ind_row_counter,1))=1;
            Aeq(Aeq_row_counter,wh_ind(wh_ind_row_counter,2))=1;
            beq(Aeq_row_counter) = 1;
            Aeq_row_counter = Aeq_row_counter + 1;
            wh_ind_row_counter=wh_ind_row_counter+1;
        end

        %>=
        for k = 1:n
            %eq1
            Ain(Ain_row_counter, T_ind(i, j)) = 1;
            Ain(Ain_row_counter, w2_ind(n*(i-1)+j, k)) = M;
            bin(Ain_row_counter) = A(i, k)+M;
            Ain_row_counter = Ain_row_counter + 1;
            %eq2
            Ain(Ain_row_counter, X_ind(k, j)) = -1;
            Ain(Ain_row_counter, T_ind(i, j)) = 1;
            Ain(Ain_row_counter, w2_ind(n*(i-1)+j, k)) = M;
            bin(Ain_row_counter) =M;
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
            %eq1
            Ain(Ain_row_counter, Y_ind(i, k)) = 1;
            Ain(Ain_row_counter, R_ind(i, j)) = -1;
            Ain(Ain_row_counter, wh_ind(wh_ind_row_counter,1))=M;
            bin(Ain_row_counter) = M;
            Ain_row_counter = Ain_row_counter + 1;
            %eq2
            Ain(Ain_row_counter, R_ind(i, j)) = -1;
            Ain(Ain_row_counter, wh_ind(wh_ind_row_counter,2))=M;
            bin(Ain_row_counter) = M-B(k, j);
            Ain_row_counter = Ain_row_counter + 1;
            %eq3
            Aeq(Aeq_row_counter,wh_ind(wh_ind_row_counter,1))=1;
            Aeq(Aeq_row_counter,wh_ind(wh_ind_row_counter,2))=1;
            beq(Aeq_row_counter) = 1;
            Aeq_row_counter = Aeq_row_counter + 1;
            wh_ind_row_counter=wh_ind_row_counter+1;
        end

        %>=
        for k = 1:n
            %eq1
            Ain(Ain_row_counter, Y_ind(i, k)) = -1;
            Ain(Ain_row_counter, R_ind(i, j)) = 1;
            Ain(Ain_row_counter, w3_ind(n*(i-1)+j, k)) = M;
            bin(Ain_row_counter) = M;
            Ain_row_counter = Ain_row_counter + 1;
            %eq2
            Ain(Ain_row_counter, R_ind(i, j)) = 1;
            Ain(Ain_row_counter, w3_ind(n*(i-1)+j, k)) = M;
            bin(Ain_row_counter) = B(k, j)+M;
            Ain_row_counter = Ain_row_counter + 1;

        end
        Aeq(Aeq_row_counter, w3_ind(n*(i-1)+j, 1:n)) = ones(1, n);
        beq(Aeq_row_counter) = 1;
        Aeq_row_counter = Aeq_row_counter + 1;
    end
end

%case4: BY=R
for i = 1:n
    for j = 1:n
        %<=
        for k = 1:n
            %eq1
            Ain(Ain_row_counter, R_ind(i, j)) = -1;
            Ain(Ain_row_counter, wh_ind(wh_ind_row_counter,1))=M;
            bin(Ain_row_counter) = -B(i, k)+M;
            Ain_row_counter = Ain_row_counter + 1;
            %eq2
            Ain(Ain_row_counter, Y_ind(k, j)) = 1;
            Ain(Ain_row_counter, R_ind(i, j)) = -1;
            Ain(Ain_row_counter, wh_ind(wh_ind_row_counter,2))=M;
            bin(Ain_row_counter) = M;
            Ain_row_counter = Ain_row_counter + 1;
            %eq3
            Aeq(Aeq_row_counter,wh_ind(wh_ind_row_counter,1))=1;
            Aeq(Aeq_row_counter,wh_ind(wh_ind_row_counter,2))=1;
            beq(Aeq_row_counter) = 1;
            Aeq_row_counter = Aeq_row_counter + 1;
            wh_ind_row_counter=wh_ind_row_counter+1;
        end

        %>=
        for k = 1:n
            %eq1
            Ain(Ain_row_counter, R_ind(i, j)) = 1;
            Ain(Ain_row_counter, w4_ind(n*(i-1)+j, k)) = M;
            bin(Ain_row_counter) = B(i, k)+M;
            Ain_row_counter = Ain_row_counter + 1;
            %eq2
            Ain(Ain_row_counter, Y_ind(k, j)) = -1;
            Ain(Ain_row_counter, R_ind(i, j)) = 1;
            Ain(Ain_row_counter, w4_ind(n*(i-1)+j, k)) = M;
            bin(Ain_row_counter) =M;
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
            %eq1
            Ain(Ain_row_counter, X_ind(i, k)) = 1;
            Ain(Ain_row_counter, wh_ind5(wh_ind5_row_counter,1))=M;
            bin(Ain_row_counter) = U(i, j)+M;
            Ain_row_counter = Ain_row_counter + 1;
            %eq2
            Ain(Ain_row_counter, wh_ind5(wh_ind5_row_counter,2))=M;
            bin(Ain_row_counter) = U(i, j)+M-W(k,l);
            Ain_row_counter = Ain_row_counter + 1;
            %eq3
            Ain(Ain_row_counter, Y_ind(l, j)) = 1;
            Ain(Ain_row_counter, wh_ind5(wh_ind5_row_counter,3))=M;
            bin(Ain_row_counter) = U(i, j)+M;
            Ain_row_counter = Ain_row_counter + 1;
            %eq4
            Aeq(Aeq_row_counter,wh_ind5(wh_ind5_row_counter,1))=1;
            Aeq(Aeq_row_counter,wh_ind5(wh_ind5_row_counter,2))=1;
            Aeq(Aeq_row_counter,wh_ind5(wh_ind5_row_counter,3))=1;
            beq(Aeq_row_counter) = 1;
            Aeq_row_counter = Aeq_row_counter + 1;
            wh_ind5_row_counter=wh_ind5_row_counter+1;
            end
        end

        %>=
        for k = 1:n
            for l=1:n
            %eq1
            Ain(Ain_row_counter, X_ind(i, k)) = -1;
            Ain(Ain_row_counter, w5_ind(n*(i-1)+j, n*(k-1)+l)) = M;
            bin(Ain_row_counter) = -U(i, j) + M;
            Ain_row_counter = Ain_row_counter + 1;
            %eq2
            Ain(Ain_row_counter, w5_ind(n*(i-1)+j, n*(k-1)+l)) = M;
            bin(Ain_row_counter) = -U(i, j) + M+W(k,l);
            Ain_row_counter = Ain_row_counter + 1;
            %eq3
            Ain(Ain_row_counter, Y_ind(l, j)) = -1;
            Ain(Ain_row_counter, w5_ind(n*(i-1)+j, n*(k-1)+l)) = M;
            bin(Ain_row_counter) = -U(i, j) + M;
            Ain_row_counter = Ain_row_counter + 1;

            end
        end
        Aeq(Aeq_row_counter, w5_ind(n*(i-1)+j, 1:n^2)) = ones(1, n^2);
        beq(Aeq_row_counter) = 1;
        Aeq_row_counter = Aeq_row_counter + 1;
    end
end

f = zeros(1, 4*n^2+12*n^3+4*n^4);
intcon = 1:4*n^2+12*n^3+4*n^4;
lb = -inf(4*n^2+12*n^3+4*n^4, 1);
%lb = zeros(4*n^2+12*n^3+4*n^4, 1)-1000;
lb(4*n^2+1:end) = 0;
ub=inf(4*n^2+12*n^3+4*n^4, 1);
%ub = zeros(4*n^2+12*n^3+4*n^4, 1)+1000;
ub(4*n^2+1:end) = 1;

model.A = sparse([Ain; Aeq]);
model.obj = f;
model.modelsense = 'min';
model.rhs = [bin; beq];
model.lb = lb;
model.ub = ub;
model.vtype = [repmat('C', 1, 4*n^2), repmat('B', 1, 12*n^3+4*n^4)];
model.sense = [repmat('<', size(Ain, 1), 1); repmat('=', size(Aeq, 1), 1)];



params.outputflag = 1;
params.BranchDir=  1;
params.MIPFocus=  0;
params.TimeLimit =60*1000;
result = gurobi(model, params);
            if strcmp(result.status, 'OPTIMAL')
                successful_trials_count = successful_trials_count + 1;
                trial_times(successful_trials_count) = result.runtime;
            else
                disp(['Trial did not find an optimal solution.']);
            end
        end
        
        % Calculate average time for successful trials
        average_times(n_index, D_index) = mean(trial_times);
    end
end


figure;
hold on;
colors = {'b', 'g', 'r', 'c', 'm', 'y'};
for n_index = 1:length(n_values)
    plot(D_values, average_times(n_index, :), '-o', 'Color', colors{n_index}, ...
        'DisplayName', ['n = ', num2str(n_values(n_index))]);
end
hold off;
xlabel('Degree');
ylabel('Average time to recover the key (seconds)');
%title('Average Time to recover the key vs Degree for Different Dimensions');
legend('show');
grid on;
