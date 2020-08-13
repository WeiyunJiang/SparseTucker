function T = smphooi_qr_p_truncatedKron4D(X,R)
    

    
    subs = X.subs;
    vals = X.vals;
    siz = X.size;

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
        Uinit{n} = ones(size(X,n),R(n))./2;
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
        intm1{4,4}{200,200} = [];
        for i = 1:4
            for j = 1:4
                intm1{i,j} = cell(siz(i),siz(j));
            end
        end
        
        for i =1:N
            order = [1 2 3 4];
            order(i) = [];
            intm2{i} = cell(siz(order(3)),siz(order(2)),siz(order(1)));
        end
        
        for mode = 1:N
            order = [1 2 3 4];
            order(mode) =[];
            
            a = order(3);
            b = order(2);
            c = order(1);
            for ind= 1:length(vals)
                coo = subs(ind,:);
                
                ind_a = coo(a);
                ind_b = coo(b);
                ind_c = coo(c);
                
                if isempty(intm1{a,b}{ind_a,ind_b})
                    intm1{a,b}{ind_a,ind_b} = kron(U{a}(ind_a,:),U{b}(ind_b,:));
                end
                if isempty(intm2{mode}{ind_a,ind_b,ind_c})
                    intm2{mode}{ind_a,ind_b,ind_c} =kron(intm1{a,b}{ind_a,ind_b}, U{c}(ind_c,:));
                end
                
                Y{mode}(coo(mode),:) = Y{mode}(coo(mode),:)+ vals(ind)*intm2{mode}{ind_a,ind_b,ind_c};
            end
            
            U{mode} = QR_pivoting(Y{mode},R(mode));
            
            
        end
       
        

        tensor_mode4 = tenmat(Y{4}, 4, [1 2 3], [R(1) R(2) R(3) siz(4)]);
        core = ttm(tensor(tensor_mode4), U{4}', 4);
        
        
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