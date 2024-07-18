function obj_val = calculate_objective(x, n, D,T)
    f = -Inf(1, n*n); 
    for ind = 1:n*n
        for i = 1:D+1
            for j = 1:D+1
                f(ind) = max(f(ind), x(i) + x(D+1 + j)+T{i,j}(ind));
            end
        end
    end
    obj_val = sum(f.^2);
end