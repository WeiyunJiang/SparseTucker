function T = smphooi_qr_p_4(X,R)
    

    
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
        Uinit{n} =ones(size(X,n),R(n))./2;
    end
    U = Uinit;
    fit = 0;
    
    if printitn > 0
        fprintf('\nTucker Kronecker Product Method qr-p:\n');
    end
    Yinitial = cell(N,1);
    for n = 1:N
        Yinitial{n} = zeros(siz(n),prod_ind(R,n));
    end
    
    

    for iter = 1:maxiters
        fitold = fit;
        Y = Yinitial;
        
        
        for mode = 1:N
            %three_ind = two_number_exc(mode);
            order = [1 2 3 4];
            order(mode) =[];
            for ind= 1:length(vals)
                coo = subs(ind,:);
                intm = kron(U{order(3)}(coo(order(3)),:),U{order(2)}(coo(order(2)),:));
                Y{mode}(coo(mode),:) = Y{mode}(coo(mode),:)+ vals(ind)*kron(intm,U{order(1)}(coo(order(1)),:));
            end
            U{mode} = QR_pivoting(Y{mode},R(mode));
            
            
        end
%         %compute TTMc for mode 1
%         
%         for ind=1:length(vals)
%             coo = subs(ind,:);
%             intm = kron(U{2}(coo(2),:),U{3}(coo(3),:));
%             Y{1}(coo(1),:) = Y{1}(coo(1),:)+ vals(ind)*kron(intm,U{4}(coo(4),:));
%         end
%         Y{1} = Y{1}*P;
%         U{1} = TRSVD(Y{1},R(1));
%         
%         %compute TTMc for mode 2
%         for ind=1:length(vals)
%             coo = subs(ind,:);
%             intm = kron(U{1}(coo(1),:),U{3}(coo(3),:));
%             Y{2}(coo(2),:) = Y{2}(coo(2),:)+ vals(ind)*kron(intm,U{4}(coo(4),:));
%         end
%         Y{2} = Y{2}*P;
%         U{2} = TRSVD(Y{2},R(2));    
%         
%         %compute TTMc for mode 3
%         for ind=1:length(vals)
%             coo = subs(ind,:);
%             intm = kron(U{1}(coo(1),:),U{2}(coo(2),:));
%             Y{3}(coo(3),:) = Y{3}(coo(3),:)+ vals(ind)*kron(intm,U{4}(coo(4),:));
%         end
%         Y{3} = Y{3}*P;
%         U{3} = TRSVD(Y{3},R(3));  
%         
%         %compute TTMc for mode 4
%         for ind=1:length(vals)
%             coo = subs(ind,:);
%             intm = kron(U{1}(coo(1),:),U{2}(coo(2),:));
%             Y{4}(coo(4),:) = Y{4}(coo(4),:)+ vals(ind)*kron(intm,U{3}(coo(3),:));
%         end
%         Y{4} = Y{4}*P;
%         U{4} = TRSVD(Y{4},R(4));        
        

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