function QRP_test
clc
A = rand(1000,1000);

tic
[Q1, R1, perm1]=qr(A,'vector');
toc

tic
[Q2,R2,perm2]=qr_Householder_pivoting_b(A);
toc
tic
[Q3,R3,perm3]=qr_Householder_pivoting_ori(A);
toc
end