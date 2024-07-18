n = 10;
mm = -1000;
mM = 1000;
D = 50;
pm = -1000;
pM = 1000;
[key,U,V,A,B,W] = GenerateKeyStickelsmaxmin(n, mm, mM, D, pm, pM);

T = cell(D + 1, D + 1);

for alpha = 0:D
    for beta = 0:D
        tem = MaxMinMulti(MaxMinpMatPower(A, alpha), W);
        T{alpha + 1, beta + 1} = MaxMinMulti(tem, MaxMinpMatPower(B, beta));
    end
end

 
objective_function = @(x) Maxmincalculate_objective(x, n, D, T,U);


lb = -1000 * ones(2 * (D + 1), 1);
ub = 1000 * ones(2 * (D + 1), 1);


options = optimoptions('simulannealbnd', 'Display', 'iter', ...
    'TemperatureFcn', @temperatureexp,'InitialTemperature',1000000,'ObjectiveLimit',0.1,'FunctionTolerance',0,'MaxStallIterations',100);

%fval_threshold = 1e-3;  
%fval = Inf;  
%while abs(fval) > fval_threshold
    % Call simulannealbnd to optimize the objective function
    x0 = mm + (mM - mm) .* rand(2 * (D + 1), 1);
    [x_opt, fval] = simulannealbnd(objective_function, x0, lb, ub, options);
    
    
    disp(['Current fval: ', num2str(fval)]);
%end

K_attack = MaxMinMulti(Applypolynomialmaxmincell([0:D; x_opt(1:D + 1)'], A), ...
    MaxMinMulti(V, Applypolynomialmaxmincell([0:D; x_opt(D + 1 + 1:end)'], B)));


key
Recovered_key=round(K_attack)
key == round(K_attack)