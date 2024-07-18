function V_X = MaxplusApplyPolynomial(V, X)
    
    [n, m] = size(X);
    [p, q] = size(V);
    D = zeros(n, m) + (-inf);
    temp = D;
    for i = 1:q
        if V(1, i) == 0
            temp = MaxplusId(n) + V(2, i);
            D = max(D, temp);
        else
            c = fastpowermaxplus(X, V(1, i));
            temp = V(2, i) + c;
            D = max(D, temp);
        end
    end
    
    V_X = D;
end