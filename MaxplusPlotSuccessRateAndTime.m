function [] = MaxplusPlotSuccessRateAndTime(max_degree, num_iterations)
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
            [key, U, V, A, B,W] = MaxplusGenerateKeyStickels(n, mm, mM, D,pm, pM);

            T = cell(D + 1, D + 1);
            for alpha = 0:D
                for beta = 0:D
                    tem = MaxplusMulti(fastpowermaxplus(A, alpha), W);
                    T{alpha + 1, beta + 1} = MaxplusMulti(tem, fastpowermaxplus(B, beta))-U;
                end
            end

            objective_function = @(x) calculate_objective(x, n, D, T);
            lb = -1000 * ones(2 * (D + 1), 1);
            ub = 1000 * ones(2 * (D + 1), 1);
            options = optimoptions('simulannealbnd', 'Display', 'iter', ...
    'TemperatureFcn', @temperatureexp,'InitialTemperature',1000000,'ObjectiveLimit',0.1,'PlotFcn',@saplotbestf);

            x0 = mm + (mM - mm) .* rand(2 * (D + 1), 1);
            tic;
            [x_opt] = simulannealbnd(objective_function, x0, lb, ub, options);

            time_consumption = toc;
            total_time = total_time + time_consumption;

            K_attack = MaxplusMulti(MaxplusApplyPolynomial([0:D; x_opt(1:D + 1)'], A), ...
                MaxplusMulti(V, MaxplusApplyPolynomial([0:D; x_opt(D + 1 + 1:end)'], B)));
            if isequal(round(K_attack), key)
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
    ylabel('Time Consumption (s)');
    title('Time Consumption for each Degree');
end
