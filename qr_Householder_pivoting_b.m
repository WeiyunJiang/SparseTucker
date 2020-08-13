function [Q,R,r]= qr_Householder_pivoting_b (A)
[m,n]=size(A);
Q=eye(m);
c=zeros(1,n);
for j =1:n
    c(j) = A(1:m,j)'*A(1:m,j);
end
r=0; tau = max(c);
[~,k]=max(c);
while tau>0
    r = r+1;
    %piv(r) = k;
    A(1:m,[r k]) = A(1:m,[k r]);
    c(:,[r k]) = c(:,[k r]);
    [v,beta]=house(A(r:m,r));
    H = (eye(m-r+1)-beta*(v*v'));
    if r <m
    A(r:m,r:n) = H*A(r:m,r:n);
    Q(:,r:m) = Q(:,r:m) -  Q(:,r:m) * (v * v' * beta);
    A(r+1:m,r) = v(2:m-r+1);
    end
    
    
    for i = r+1:n
        c(i) = c(i)-A(r,i)^2;
    end
    if r<n
        tau = max(c(r+1:n));
        [~,k]=max(c(r+1:n));
        k = k+r;
    else
        tau = 0;
    end
end
R = triu(A);

end