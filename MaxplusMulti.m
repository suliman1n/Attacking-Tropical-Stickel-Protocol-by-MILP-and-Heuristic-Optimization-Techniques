function C = MaxplusMulti(A,B)

[m,n]=size(A);
[k,l]=size(B);
if (n~=k)  
    error("Dimension Error!")
end
     Ax = [];     
    Ax1 = []; 
for i=1:m
    for j=1:l
        Ax1 = A(i,:) + B(:,j)';
        Ax(i,j) = max(Ax1);
        
    end
end  
C=Ax;
