function phi = THOMAS(main,sup,sub,rhs)

% main = main diag
% sup = super (above) diag
% sub = sub (below) diag
% rhs = source (RHS) vector 

% Pass in entire vectors from coefficient matrix (including BC points)
n=length(rhs);
phi = zeros(length(rhs),1);

%Forward Elimination Loop
for i=2:n
 xmlt=sub(i-1)/main(i-1);
 main(i) = main(i) - xmlt*sup(i-1);
 rhs(i) = rhs(i) - xmlt*rhs(i-1);
end

% Back Substitution
phi(n) = rhs(n)/main(n);
for i=n-1:-1:1
 phi(i) = (rhs(i) - sup(i)*phi(i+1))/main(i);
end
end