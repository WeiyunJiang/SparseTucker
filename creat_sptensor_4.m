function X =  creat_sptensor_4 (s)
R=10;
A = sprand(700,R,s);
B = sprand(800,R,s);
C = sprand(500,R,s);
D = sprand(20,R,s);
%G = tenrand([100 100 100 100]);
G =sptenrand([R R R R], round(R^4*0.001));
T = ttm(G, A, 1); 
T = ttm(T, B, 2); 
T = ttm(T, C, 3); 
T = ttm(T, D, 4); 
S = find(T > 1e-4);
V = T(S);
X = sptensor(S,V,[700 800 500 20]);
sparsity = nnz(X)/(prod(size(X)));
disp(['sparsity is ',num2str(sparsity)]);
end