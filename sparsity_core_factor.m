function [sparsity_core,sparsity_fac]=sparsity_core_factor(A,R)
nnz = length(find(A.core>0.00005 | A.core<-0.00005 ));
sparsity_core = nnz/prod(size(A.core));
for i = 1:R
    nnz = length(find(A.U{R}>0.00005 | A.U{R}<-0.00005 ));
    sparsity_fac{i} = nnz/prod(size(A.U{i}));
end
end