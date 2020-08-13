function U=  QR_pivoting(Y, R_n)
[Q,R,P] = qr(Y*Y');
%[Q,R,r]= qr_Householder_pivoting_b (Y*Y');
U = Q(:,1:R_n);
end

