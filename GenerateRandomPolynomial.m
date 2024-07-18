function p = GenerateRandomPolynomial(D, pm, pM)
    %d = randi([1, D]);
    d=D;
    p = zeros(2, d + 1); 

    for i = 1:(d + 1)
        p(:, i) = [i - 1; randi([pm, pM])]; 
    end
end