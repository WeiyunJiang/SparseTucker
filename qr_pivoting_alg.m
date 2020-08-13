function [Q,R,P] = qr_pivoting_alg(A)
[m,n]=size(A);
R = A;
Q = eye(m);

% I am going to use a permutation matrix.
P = eye(n);

% Compute the norms.
for i = 1 : n
  colnorm(i) = R(:,i)'*R(:,i)
end

%Swapping procedure
for i=1:n-1

   %Find max col norm
   max_colnorm = colnorm(i); perms = i;
   for  j = i + 1 : n
       if (colnorm(j) > max_colnorm)
           perms = j;
           max_colnorm = colnorm(j);
       end
   end

   %Break 
   if ( colnorm(perms) == 0 )
       break;
   end

   %Swap P
   temp = P(:, i);
   P(:, i) = P (:, perms)
   P(:, perms) = temp

   %Swap R
   temp = R(:, i);
   R(:, i) = R (:, perms)
   R(:, perms) = temp

   %Swap colnorm
   colnorm = colnorm*P

   % Get the Householder vector from get_house.
   v = get_house(R(:,i),i,m)

   % Apply the transformation to R from the left.
   R = R-v*(v'*R)

   % And also apply it to Q from the right.
   Q = Q - (Q*v)*v'

   %Norm downdate
   if i~=(n-1)
       colnorm(i+1:n) = colnorm(i+1:n) - R(i, i+1 : n).^2 
   end

end

% Get the Householder vector from get_house.
v = get_house(R(:,n),n,m)

% Apply the transformation to R from the left.
R = R- v  *(v' * R)

% And also apply it to Q from the right.
Q = Q - (Q * v)* v'
end