function T = smphooi_qr_p_parfor_truncatedKron3D(X,R)
    

    
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
        fprintf('\nTucker Kronecker Product Method qr-p parfor truncated:\n');
    end
    Yinitial = cell(N,1);
    for n = 1:N
        Yinitial{n} = zeros(siz(n),prod_ind(R,n));
    end
    [ul,J]=symbolicTTMC(subs,vals,N);
  

    for iter = 1:maxiters
        fitold = fit;
        
        for mode = 1:N
            Y = Yinitial;
            order = [1 2 3];
            order(mode) =[]; 
            temp = zeros(length(J{mode}),prod_ind(R,mode));
            parfor i_ind= 1:length(J{mode})
                %i = zeros(length(J{mode}),1);
                i(i_ind,1) = J{mode}{i_ind}(1,1);
                
                temp(i_ind,:) = zeros(1,prod_ind(R,mode));
                
                for nz = 1:length(ul{mode}{i_ind})
                    row = J{mode}{i_ind}(nz,:);
                    
                    temp(i_ind,:) = temp(i_ind,:)+ ul{mode}{i_ind}(nz)*kron(U{order(2)}(row(3),:),U{order(1)}(row(2),:));
                end
            end
            
            Y{mode}(i(1:length(J{mode})),:) = temp(1:length(J{mode}),:);
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