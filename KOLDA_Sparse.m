X = rand(5,6,4,2); R = [2 3]; C = [4 1];
I = size(X); J = prod(I(R)); K = prod(I(C));
Y = reshape(permute(X,[R C]),J,K); % convert X to matrix Y
Z = ipermute(reshape(Y,[I(R) I(C)]),[R C]); % convert back to tensor