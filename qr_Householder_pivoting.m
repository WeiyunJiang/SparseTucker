% Variant of Businger and Golub's QR with column pivoting
% Copyright (c) 2016 by Pranay Seshadri
% norm is positive
function [Q, R,p] = qr_Householder_pivoting(A)
R=A;
[m,n] = size(A);
p = 1 : n;
lens = zeros(n,1); % Initialize column norms vector
Q = eye(m,m);

% Compute column norms
for j = 1 : n
    lens(j) = R(:,j)'*R(:,j);
end

% Reduction steps
for k = 1 : min(m,n) - 1
    
    % Compute max column norm
    [~,j_star] = max(lens(p(k:n)));
    j_star = j_star + (k - 1);
    
    % Swap jth and pth columns!
    if(k ~= j_star)
        
%         temp = R(1:m,k);
%         R(1:m,k) = R(1:m,j_star);
%         R(1:m,j_star) = temp;
%         clear temp;
        
        temp = p(k);
        p(k) = p(j_star);
        p(j_star) = temp;
        clear temp;
        
    end
    v = get_house(R(:,p(k)),k,m);
    R = R-v*(v'*R);
    Q = Q - (Q*v)*v';
    % Reduction -- compute Householder matrix
%     [v,betav] = house(A(k:m,k));
%     H = (eye(m-k+1) - betav * (v * v') ); % I'd prefer not to use H{j}!
%     A(k:m,k:n) =  H * A(k:m,k:n);
%     Q(:,k:m) = Q(:,k:m) -  Q(:,k:m) * (v * v' * betav);
    
    % Update the remaining column norms
    if(k~=n)
        for j = k + 1 : n
            lens(p(j)) = lens(p(j)) - R(1,p(j))^2;
        end
    end
end
v = get_house(R(:,p(min(m,n))),3,6);
R = R-v*(v'*R);
Q = Q - (Q*v)*v';
R = R(:,p);


end