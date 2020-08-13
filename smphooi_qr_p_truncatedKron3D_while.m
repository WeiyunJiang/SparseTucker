function T = smphooi_qr_p_truncatedKron3D_while(X,R)
    

    
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
        Uinit{n} = rand(size(X,n),R(n));
    end
    U = Uinit;
    fit = 0;
    
    if printitn > 0
        fprintf('\nTucker Kronecker Product Method qr-p truncated while:\n');
    end
    Yinitial = cell(N,1);
    for n = 1:N
        Yinitial{n} = zeros(siz(n),prod_ind(R,n));
    end
    
    
    
    for iter = 1:maxiters
        fitold = fit;
        Y = Yinitial;
        
        
        for mode = 1:N
            U_ind=two_number_exc_3(mode);
            nnz = length(vals);
            subs_0 = subs;
            vals_0 = vals;
            col_mode = subs(:,mode);
            subs_0 (:,mode) = [];
            subs_0 = [col_mode subs_0];
            while nnz>0
                index = subs_0(1,:);
                i_2=subs_0(:,2);
                i_3=subs_0(:,3);
                index_share_2 = find(i_2 == index(2));
                z = i_3(index_share_2)==index(3);
                index_share = index_share_2(z);
                
                intm = kron(U{U_ind(2)}(index(3),:),U{U_ind(1)}(index(2),:));
                
                Y{mode}(subs_0(index_share,1),:) = Y{mode}(subs_0(index_share,1),:)+ vals_0(index_share)*intm;
                
                subs_0(index_share,:)= [];
                vals_0(index_share) = [];
                nnz = nnz - length(index_share);
                
            end
            
            U{mode} = QR_pivoting(Y{mode},R(mode));
            
            
        end
        
        

        tensor_mode3 = tenmat(Y{3}, 3, [1 2 ], [R(1) R(2) siz(3)]);
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