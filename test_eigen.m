Y = randn(5,5);
R=5;
U1 = TRSVD(Y, R);
[U2,D] = eig(Y);
U3 =  QR_pivoting(Y, R);