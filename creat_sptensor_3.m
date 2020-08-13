function Y =  creat_sptensor_3 (s)
R=50;
A = sprand(200,R,s);
B = sprand(200,R,s);
C = sprand(200,R,s);

%G = tenrand([100 100 100 100]);
G =sptenrand([R R R], round(R^3*0.01));
T = ttm(G, A, 1); 
T = ttm(T, B, 2); 
T = ttm(T, C, 3); 

S = find(T > 1e-4);
V = T(S);
Y = sptensor(S,V,[200 200 200]);
sparsity = nnz(Y)/(prod(size(Y)));
disp(['sparsity is ',num2str(sparsity)]);
end