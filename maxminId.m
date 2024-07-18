function A = maxminId(n)
%TropIdent - Creates a Tropical Identity matrix of size n 
%   Detailed explanation goes here
A = zeros(n);
A(1:n,1:n) = -inf; %-inf
for i=1:n
    A(i,i)=inf; %inf
end


%min(B(A > B))
%find(A >= B & c >= B)