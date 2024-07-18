n = 5;
mm = -1000;
mM = 1000;
D = 5;
pm = -1000;
pM = 1000;
[key, U, V, A, B,W] = MaxplusGenerateKeyStickels(n, mm, mM, D, pm, pM);

T = cell(D + 1, D + 1);

for alpha = 0:D
    for beta = 0:D
        tem = MaxplusMulti(fastpowermaxplus(A, alpha), W);
        T{alpha + 1, beta + 1} = MaxplusMulti(tem, fastpowermaxplus(B, beta))-U;
    end
end



Ain = [];
            bin = [];
            Aeq = [];
            beq = [];
            M = 1000000;
            %< side
            for ind = 1:n*n
                for i = 1:D+1
                    for j = 1:D+1
                        Ain = [Ain; zeros(1, D+1+D+1+((D+1)^2)*n^2)];
                        Ain(end, i) = 1;
                        Ain(end, D+1+j) = 1;
                        bin = [bin; -T{i,j}(ind)];
                    end
                end
            end

            %> side
            for ind = 1:n*n
                w_counter = 1;
                for i = 1:D+1
                    for j = 1:D+1
                        Ain = [Ain; zeros(1, D+1+D+1+((D+1)^2)*n^2)];
                        Ain(end, i) = -1;
                        Ain(end, D+1+j) = -1;
                        Ain(end, D+1+D+1+w_counter+((D+1)^2)*(ind-1)) = M;
                        bin = [bin; T{i,j}(ind) + M];
                        w_counter = w_counter + 1;
                    end
                end
                Aeq = [Aeq; zeros(1, D+1+D+1+((D+1)^2)*n^2)];
                Aeq(end, D+1+D+1+1+((D+1)^2)*(ind-1):D+1+D+1+((D+1)^2)+((D+1)^2)*(ind-1)) = ones(1, (D+1)^2);
                beq = [beq; 1];
            end

            f = zeros(1, D+1+D+1+((D+1)^2)*n^2);
            intcon = 1:D+1+D+1+((D+1)^2)*n^2;
            lb = zeros(D+1+D+1+((D+1)^2)*n^2, 1) - inf;
            lb(D+1+D+1+1:end) = 0;
            ub = zeros(D+1+D+1+((D+1)^2)*n^2, 1) + inf;
            ub(D+1+D+1+1:end) = 1;

            model.A = sparse([Ain; Aeq]);
model.rhs = [bin; beq];
model.sense = [repmat('<', size(Ain, 1), 1); repmat('=', size(Aeq, 1), 1)];
model.obj = f;
model.modelsense = 'min';
model.lb = lb;
model.ub = ub;
model.vtype = [repmat('C', 1, D+1+D+1), repmat('B', 1, ((D+1)^2)*n^2)];


params.outputflag = 1;
params.BranchDir=  -1;
params.MIPFocus=  1;


result = gurobi(model, params);

    if strcmp(result.status, 'OPTIMAL')
    z = result.x;
    x = z(1:D+1);
    y = z(D+1+1:D+1+D+1);
else
    x = [];
    y = [];
end



    K_attack = MaxplusMulti(MaxplusApplyPolynomial([0:D;x'], A),MaxplusMulti(V, MaxplusApplyPolynomial([0:D; y'], B)))
    key==K_attack
    