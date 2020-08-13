function U = LMSVD_R(Y, R)
    Z = Y*Y';
    r= R;
    opts.tol = 1e-8;
    opts.maxit = 150;
    [U,S,V] = lmsvd(Z,r, opts);
    %[eigvec,piv] = sort(diag(S),'descend');
    
    %U = U1(:,piv(1:R));
end