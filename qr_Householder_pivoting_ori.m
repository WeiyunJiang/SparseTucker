% Variant of Businger and Golub's QR with column pivoting
% Copyright (c) 2016 by Pranay Seshadri
function [Q, R,perm] = qr_Householder_pivoting_ori(A)
[m,n] = size(A);
perm = 1 : n;
c = zeros(1,n); % Initialize column norms vector
Q = eye(m,m);

% Compute column norms
for j =1:n
    c(j) = A(1:m,j)'*A(1:m,j);
end

% Reduction steps
for r = 1 : min(m,n) - 1
    
    % Compute max column norm
    [~,j_star] = max(c(r:n));
    j_star = j_star + (r - 1);
    
    % Swap jth and pth columns!
    if(r ~= j_star)
        
        temp = A(1:m,r);
        A(1:m,r) = A(1:m,j_star);
        A(1:m,j_star) = temp;
        clear temp;
        
        temp = perm(r);
        perm(r) = perm(j_star);
        perm(j_star) = temp;
        clear temp;
        
    end
    
    % Reduction -- compute Householder matrix
    [v,beta] = house(A(r:m,r));
    H = (eye(m-r+1) - beta * (v * v') ); % I'd prefer not to use H{j}!
    A(r:m,r:n) =  H * A(r:m,r:n);
    Q(:,r:m) = Q(:,r:m) -  Q(:,r:m) * (v * v' * beta);
    A(r+1:m,r) = v(2:m-r+1);
    % Update the remaining column norms
    if(r~=n)
        for i = r+1:n
            c(i) = c(i)-A(r,i)^2;
        end
    end
end


R = triu(A);

end