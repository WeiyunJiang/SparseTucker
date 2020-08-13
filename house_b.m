function [v,beta] = house_b(x)
n = length(x);
sigma = x(2:n)'*x(2:n);
v = [1;x(2:n)];
if sigma == 0
    beta = 0;
else
    mu = sqrt(x(1)^2+sigma);
    if x(1) <= 0
        v(1)=x(1)-mu;
    else
        v(1) = -sigma/(x(1)+mu);
    end
    beta =2*v(1)^2/(sigma+v(1)^2);
    v=v/v(1);
end