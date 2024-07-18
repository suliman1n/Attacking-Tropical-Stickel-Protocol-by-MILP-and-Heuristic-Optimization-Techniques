function [C] = MaxMinpMatPower(A,p)

[n,m] = size(A);
if (n~=m)
    error("Dimension Error! Not a square matrix")
end 
if p ==0
    C = maxminId(n);
    return
end 
temp = A;
for i =2:p
    temp = MaxMinMulti(A,temp);
end
C = temp;
end