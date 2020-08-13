function U = SSVD(Y, R)
    Z = Y*Y';
    [U1,S,V] = svds(Z);
    [eigvec,piv] = sort(diag(S),'descend');

    U = -U1(:,piv(1:R));
end