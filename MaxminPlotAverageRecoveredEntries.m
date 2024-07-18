function [] = MaxminPlotAverageRecoveredEntries(max_degree, num_iterations)
    avg_recovered_entries = zeros(1, max_degree - 1);
    time_consumptions = zeros(1, max_degree - 1);

    for D = 2:max_degree
        total_recovered_entries = 0;
        total_time = 0;
        
        for i = 1:num_iterations
            n = 10;
            mm = -1000;
            mM = 1000;
            pm = -1000;
            pM = 1000;
            [key, U, V, A, B,W] = GenerateKeyStickelsmaxmin(n, mm, mM, D, pm, pM);

            T = cell(D + 1, D + 1);

            for alpha = 0:D
                for beta = 0:D
                    tem = MaxMinMulti(MaxMinpMatPower(A, alpha), W);
                    T{alpha + 1, beta + 1} = MaxMinMulti(tem, MaxMinpMatPower(B, beta));
                end
            end

            objective_function = @(x) Maxmincalculate_objective(x, n, D, T, U);

            lb = -1000 * ones(2 * (D + 1), 1);
            ub = 1000 * ones(2 * (D + 1), 1);

            options = optimoptions('simulannealbnd', 'Display', 'iter', ...
                'TemperatureFcn', @temperatureexp, 'InitialTemperature', 1000000, 'ObjectiveLimit', 0.1, 'FunctionTolerance', 0, 'MaxStallIterations', 300);

            x0 = mm + (mM - mm) .* rand(2 * (D + 1), 1);
            tic;
            [x_opt] = simulannealbnd(objective_function, x0, lb, ub, options);

            time_consumption = toc;
            total_time = total_time + time_consumption;

            K_attack = MaxMinMulti(Applypolynomialmaxmincell([0:D; x_opt(1:D + 1)'], A), ...
                MaxMinMulti(V, Applypolynomialmaxmincell([0:D; x_opt(D + 1 + 1:end)'], B)));
            matches = sum(round(K_attack) == key, 'all');
            total_recovered_entries = total_recovered_entries + matches;
        end
        avg_recovered_entries(D - 1) = total_recovered_entries / num_iterations;
        time_consumptions(D - 1) = total_time / num_iterations;
    end

    
    figure;
    plot(2:max_degree, avg_recovered_entries, '-o', 'LineWidth', 2);
    xlabel('Degree');
    ylabel('Average Recovered Entries');
    title('Average Number of Recovered Entries for each Degree');

    
    figure;
    plot(2:max_degree, time_consumptions, '-o', 'LineWidth', 2);
    xlabel('Degree');
    ylabel('Time Consumption (Seconds)');
    title('Time Consumption for each Degree');
end
