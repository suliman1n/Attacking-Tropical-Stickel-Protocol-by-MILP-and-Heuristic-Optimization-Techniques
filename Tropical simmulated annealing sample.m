n = 10;
mm = -1000;
mM = 1000;
D = 10;
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

%x0 = zeros(2 * (D + 1), 1);  % Initial guess for the variables x 
objective_function = @(x) calculate_objective(x, n, D, T);


lb = -1000 * ones(2 * (D + 1), 1); % Lower bounds
ub = 1000 * ones(2 * (D + 1), 1); % Upper bounds


options = optimoptions('simulannealbnd', 'Display', 'iter', ...
    'TemperatureFcn', @temperatureexp,'InitialTemperature',1000000,'ObjectiveLimit',0.1,'PlotFcn',@saplotbestf);


%fval_threshold = 1e-3;  
%fval = Inf;  
%while abs(fval) > fval_threshold
    x0 = mm + (mM - mm) .* rand(2 * (D + 1), 1);
    [x_opt, fval] = simulannealbnd(objective_function, x0, lb, ub, options);
    
    
    disp(['Current fval: ', num2str(fval)]);
%end

K_attack = MaxplusMulti(MaxplusApplyPolynomial([0:D; x_opt(1:D + 1)'], A), ...
    MaxplusMulti(V, MaxplusApplyPolynomial([0:D; x_opt(D + 1 + 1:end)'], B)));

key
Recovered_key=round(K_attack)
%key == round(K_attack)
%fval