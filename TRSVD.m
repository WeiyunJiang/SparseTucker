function U = TRSVD(Y, R)
    Z = Y*Y';
    [V,D] = eig(Z);
    [eigvec,pi] = sort(diag(D),'descend');

    U = V(:,pi(1:R));
end