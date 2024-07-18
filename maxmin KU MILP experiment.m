max_degree=10;
num_iterations=1;
    success_rates = zeros(1, max_degree - 1);
    time_consumptions = zeros(1, max_degree - 1);

    for D = 2:max_degree
        successful_attempts = 0;
        total_time = 0;

        for i = 1:num_iterations
            n = 10;
            mm = -1000;
            mM = 1000;
            pm = -1000;
            pM = 1000;
            [key,U,V,A,B,W] = GenerateKeyStickelsmaxmin(n, mm, mM, D, pm,pM);
            tic;
            w_h_ind=reshape(2*(D+1)+n^2*(D+1)^2+1:2*(D+1)+n^2*(D+1)^2+3*n^2*(D+1)^2, 3, []).';
T = cell(D + 1, D + 1);

for alpha = 0:D
    for beta = 0:D
        tem = MaxMinMulti(MaxMinpMatPower(A, alpha), W);
        T{alpha + 1, beta + 1} = MaxMinMulti(tem, MaxMinpMatPower(B, beta));
    end
end



Ain=[];
bin=[];
Aeq=[];
beq=[];
M=10000;
w_h_row_counter=1;
%< side
for ind = 1:n*n
    for i = 1:D+1
        for j = 1:D+1
            Ain=[Ain;zeros(3,D+1+D+1+((D+1)^2)*n^2+3*n^2*(D+1)^2)];
            %eq1
            Ain(end-2,i)=1;
            Ain(end-2,w_h_ind(w_h_row_counter,1))=M;
            bin = [bin; U(ind)+M];
            %eq2
            Ain(end-1,D+1+j)=1;
            Ain(end-1,w_h_ind(w_h_row_counter,2))=M;
            bin = [bin; U(ind)+M];
            %eq3
        
            Ain(end,w_h_ind(w_h_row_counter,3))=M;
            bin = [bin; M-T{i,j}(ind)+U(ind)];
            %eq4
            Aeq=[Aeq;zeros(1,D+1+D+1+((D+1)^2)*n^2+3*n^2*(D+1)^2)];
            Aeq(end,w_h_ind(w_h_row_counter,1))=1;
            Aeq(end,w_h_ind(w_h_row_counter,2))=1;
            Aeq(end,w_h_ind(w_h_row_counter,3))=1;
            beq = [beq; 1];
            w_h_row_counter=w_h_row_counter+1;
        end
    end
end


%> side
for ind = 1:n*n
    w_counter=1;
    for i = 1:D+1
        for j = 1:D+1
            Ain=[Ain;zeros(3,D+1+D+1+((D+1)^2)*n^2+3*n^2*(D+1)^2)];
            %eq1
            Ain(end-2,i)=-1;
            Ain(end-2,D+1+D+1+w_counter+((D+1)^2)*(ind-1))=M;
            bin = [bin; -U(ind)+M];
            %eq2
            Ain(end-1,D+1+j)=-1;
            Ain(end-1,D+1+D+1+w_counter+((D+1)^2)*(ind-1))=M;
            bin = [bin; -U(ind)+M];
            %eq3
            
            Ain(end,D+1+D+1+w_counter+((D+1)^2)*(ind-1))=M;
            bin = [bin; -U(ind)+M+T{i,j}(ind)];

            w_counter=w_counter+1;
        end
    end
    Aeq=[Aeq;zeros(1,D+1+D+1+((D+1)^2)*n^2+3*n^2*(D+1)^2)];
    Aeq(end,D+1+D+1+1+((D+1)^2)*(ind-1):D+1+D+1+((D+1)^2)+((D+1)^2)*(ind-1))=ones(1, (D+1)^2);
    beq=[beq;1];

end



    f = zeros(1, D+1+D+1+((D+1)^2)*n^2+3*n^2*(D+1)^2);

    intcon = 1:D+1+D+1+((D+1)^2)*n^2+3*n^2*(D+1)^2;

    lb=zeros(D+1+D+1+((D+1)^2)*n^2+3*n^2*(D+1)^2,1)-inf;
    lb(D+1+D+1+1:end)=0;

    ub=zeros(D+1+D+1+((D+1)^2)*n^2+3*n^2*(D+1)^2,1)+inf;
    ub(D+1+D+1+1:end)=1;

    
    model.A = sparse([Ain; Aeq]);
model.rhs = [bin; beq];
model.sense = [repmat('<', size(Ain, 1), 1); repmat('=', size(Aeq, 1), 1)];
model.obj = f;
model.modelsense = 'min';
model.lb = lb;
model.ub = ub;
model.vtype = [repmat('C', 1, D+1+D+1), repmat('B', 1, ((D+1)^2)*n^2+3*n^2*(D+1)^2)];

params.outputflag = 1;
params.BranchDir=  1;
params.MIPFocus=  0;
params.TimeLimit =60*1000;
result = gurobi(model, params);

if strcmp(result.status, 'OPTIMAL')
    z = result.x;
    x = z(1:D+1);
    y = z(D+1+1:D+1+D+1);
else
    x = []
    y = []
end



            K_attack = MaxMinMulti(Applypolynomialmaxmincell([0:D;x'], A),MaxMinMulti(V, Applypolynomialmaxmincell([0:D; y'], B)));
            time_consumption = toc;
            total_time = total_time + time_consumption;
            if isequal(K_attack, key)
                successful_attempts = successful_attempts + 1;
            end
        end
        success_rates(D - 1) = successful_attempts / num_iterations;
        time_consumptions(D - 1) = total_time / num_iterations;
    end

    
    figure;
    plot(2:max_degree, success_rates, '-o', 'LineWidth', 2);
    xlabel('Degree');
    ylabel('Success Rate');
    title('Success Rate for each Degree');

    
    figure;
    plot(2:max_degree, time_consumptions, '-o', 'LineWidth', 2);
    xlabel('Degree');
    ylabel('Time Consumption (seconds)');
    title('Time Consumption for each Degree');

