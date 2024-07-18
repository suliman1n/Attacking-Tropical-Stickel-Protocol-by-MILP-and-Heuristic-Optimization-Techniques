function [key,U,V,A,B,W] = MaxplusGenerateKeyStickels(n, mm, mM, D, pm, pM)
    
    A = randi([mm, mM], n);
    B = randi([mm, mM], n);
    W=randi([mm, mM], n);

    
    p1 = GenerateRandomPolynomial(D, pm, pM);
    p2 = GenerateRandomPolynomial(D, pm, pM);
    q1 = GenerateRandomPolynomial(D, pm, pM);
    q2 = GenerateRandomPolynomial(D, pm, pM);

    
    temp=MaxplusMulti(MaxplusApplyPolynomial(p1, A), W);
    U=MaxplusMulti(temp, MaxplusApplyPolynomial(p2, B));

    temp=MaxplusMulti(MaxplusApplyPolynomial(q1, A), W);
    V=MaxplusMulti(temp, MaxplusApplyPolynomial(q2, B));
    
    

    
    KA = MaxplusMulti(MaxplusApplyPolynomial(p1, A), MaxplusMulti(V, MaxplusApplyPolynomial(p2, B)));
    KB = MaxplusMulti(MaxplusApplyPolynomial(q1, A), MaxplusMulti(U, MaxplusApplyPolynomial(q2, B)));

    
    if isequal(KA, KB)
        key = KA;
    else
        key = [];
    end
end
