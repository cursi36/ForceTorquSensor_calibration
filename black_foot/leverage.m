function h = leverage(V)

[r,c] = size(V);
X = zeros(r,c+1);

X(:,1) = ones(r,1);
X(:,2:1+c) = V;

H = X*(X'*X)^-1*X';
h = diag(H);




end