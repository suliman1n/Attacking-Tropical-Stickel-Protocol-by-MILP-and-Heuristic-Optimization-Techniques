function [key,U,V,A,B,W] = GenerateKeyStickelsmaxmin(n, mm, mM, D, pm, pM)
    
    A = randi([mm, mM], n);
    B = randi([mm, mM], n);
    W=randi([mm, mM], n);

    
    p1 = GenerateRandomPolynomial(D, pm, pM);
    p2 = GenerateRandomPolynomial(D, pm, pM);
    q1 = GenerateRandomPolynomial(D, pm, pM);
    q2 = GenerateRandomPolynomial(D, pm, pM);

    
    

    tem = MaxMinMulti(Applypolynomialmaxmincell(p1, A), W);
    U = MaxMinMulti(tem, Applypolynomialmaxmincell(p2, B));
    
    tem = MaxMinMulti(Applypolynomialmaxmincell(q1, A), W);
    V = MaxMinMulti(tem, Applypolynomialmaxmincell(q2, B));
    
    KA = MaxMinMulti(Applypolynomialmaxmincell(p1, A), MaxMinMulti(V, Applypolynomialmaxmincell(p2, B)));
    KB = MaxMinMulti(Applypolynomialmaxmincell(q1, A), MaxMinMulti(U, Applypolynomialmaxmincell(q2, B)));

    
    if isequal(KA, KB)
        key = KA;
    else
        key = [];
    end
end
