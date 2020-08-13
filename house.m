function [v,betav] = house(x)
n = length(x);
sigmav = x(2:n)' * x(2:n);
v = [1; x(2:n) ];

if(length(x) == 1)
    v = 1;
    betav = x(1);
elseif(sigmav == 0 && x(1) >= 0 )
    betav = 0;
elseif(sigmav == 0 && x(1) < 0 )
    betav = -2;
else
    muv = sqrt(x(1)^2 + sigmav);
    if x(1) <= 0
        v(1) = x(1) - muv;
    else
        v(1) = -sigmav/(x(1) + muv) ;
    end
    betav = 2*v(1)^2 / (sigmav + v(1)^2);
    v = v/v(1);
end