function T = smphooi_qr_p_truncatedKron3D(X,R)
    

    
    subs = X.subs;
    vals = X.vals;
    siz=X.size;

    N = ndims(X);
    normX = norm(X);
    maxiters = 50;
    fitchangetol=1e-4;
    printitn =1;
    if numel(R) == 1
        R = R * ones(N,1);
    end
    
    U = cell(N,1);
    dimorder = 1:N;
    Uinit = cell(N,1);
    for n = dimorder(2:end)
        Uinit{n} = zeros(size(X,n),R(n));
    end
    U = Uinit;
    fit = 0;
    
    if printitn > 0
        fprintf('\nTucker Kronecker Product Method qr-p truncated:\n');
    end
    Yinitial = cell(N,1);
    for n = 1:N
        Yinitial{n} = zeros(siz(n),prod_ind(R,n));
    end
    


    for iter = 1:maxiters
        fitold = fit;
        Y = Yinitial;
        intm = cell(N,1);
        for n = 1:N
            order = [1 2 3];
            order(n) =[];
            intm{n} = cell(siz(order(2)),siz(order(1)));
        end
        
        for mode = 1:N
            three_ind = two_number_exc_3(mode);
%             order = [1 2 3];
%             order(mode) =[];
            for ind= 1:length(vals)
                coo = subs(ind,:);
                ind_a = coo(three_ind(2));
                ind_b = coo(three_ind(1));   
                if isempty(intm{mode}{ind_a,ind_b})
                intm{mode}{ind_a,ind_b} = kron(U{three_ind(2)}(ind_a,:),U{three_ind(1)}(ind_b,:));
                Y{mode}(coo(mode),:) = Y{mode}(coo(mode),:)+ vals(ind)*intm{mode}{ind_a,ind_b};
                else 
                Y{mode}(coo(mode),:) = Y{mode}(coo(mode),:)+ vals(ind)*intm{mode}{ind_a,ind_b};
                end
            end
            U{mode} = QR_pivoting(Y{mode},R(mode));  
        end


        tensor_mode3 = tenmat(Y{3}, 3, [1 2], [R(1) R(2) siz(3)]);
        core = ttm(tensor(tensor_mode3), U{3}', 3);
        
        
        normresidual = sqrt( normX^2 - norm(core)^2 );
        fit = 1 - (normresidual / normX); %fraction explained by model
        fitchange = abs(fitold - fit);
        if mod(iter,printitn)==0
            fprintf(' Iter %2d: fit = %e fitdelta = %7.1e\n', iter, fit, fitchange);
        end
        if (iter > 1) && (fitchange < fitchangetol)
            break;
        end
    end
    T = ttensor(core, U);
end