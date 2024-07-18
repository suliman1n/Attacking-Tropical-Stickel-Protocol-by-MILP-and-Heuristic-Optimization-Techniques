function A = fastpowermaxplus(B, t)
    [n, m] = size(B);
    
    if (n ~= m)
        error("Dimension Error! Not a square matrix")
    end
    
    if t == 0
        A = MaxplusId(n);
        return
    end 
    
    exp = dec2bin(t);
    value = MaxplusId(n); 
    
    for i = 1:length(exp)
        value = MaxplusMulti(value, value);
        
        if exp(i) == '1'
            value = MaxplusMulti(value, B);
        end
    end
    
    A = value;
end
