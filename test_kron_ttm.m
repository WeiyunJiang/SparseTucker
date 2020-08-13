function test_kron_ttm
R= [2 2 2 2];

  subs = [1 1 1 2; 1 1 2 2; 2 2 2 1; 1 1 1 1; 1 2 1 2; 1 1 1 2];
   vals = [5; 2; 3; 4; 5; 6];
  siz = [2 2 2 2];
    X = sptensor(subs,vals,siz);
U1 = reshape (1:4,2,2);
U2 = reshape (1:4,2,2);
U3 = reshape(1:4,2,2);
U4 = reshape(1:4,2,2);
U1T = U1';
U2T=U2';
U3T=U3';
U4T = U4';
X1 = ttm(X,U1T,1);
X3 = ttm(X1,U3T,3);
X4 = ttm(X3,U4T,4);
X4_mat = tenmat(X4,2);



Y = zeros(2,8);
for ind=1:length(vals)
    coo = subs(ind,:);
    intm = kron(U4(coo(4),:),U3(coo(3),:));
    Y(coo(2),:) = Y(coo(2),:)+ vals(ind)*kron(intm,U1(coo(1),:));
end

end