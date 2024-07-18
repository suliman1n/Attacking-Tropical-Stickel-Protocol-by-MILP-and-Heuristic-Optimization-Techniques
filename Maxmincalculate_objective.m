function obj_val = Maxmincalculate_objective(x, n, D,T,U)
    f = -ones(1, n*n)*1000; %-inf
    for ind = 1:n*n
        for i = 1:D+1
            for j = 1:D+1
                f(ind) = max(f(ind), min([x(i),x(D+1 + j),T{i,j}(ind)]))-U(ind);
            end
        end
    end
    obj_val = sum(f.^2);
end