function U = TTMc(n,Jn,uln,R, Uinit,vnz)
% Jn is the values of the nonzeros
% uln is the list of the nonzeros




parfor i = 1:size(Jn)
Y = zeros(R(1),prod(R(2:end)));
    for j = 1:uln(i)
        Y(i,:) = Y(i,:) + 
        
    end
end

U{n} = TRSVD(Y,R(n));
end