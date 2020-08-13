function [ul,J]=symbolicTTMC(subs,vals,N)
    parfor n = 1:N
        order = [1 2 3];
        order(n) =[]; 
        subs_0 = subs(:,[n order]);
        vals_0 = vals;
        num=size(subs_0,1);
        i_ind = 1;
        
        while num > 0
            i_0 = subs_0(1,1);
            ind_share= find(subs_0(:,1)==i_0);
            J{n}{i_ind} = subs_0(ind_share,:);
            ul{n}{i_ind} = vals_0(ind_share);
            i_ind=i_ind+1;
            subs_0(ind_share,:) = [];
            vals_0(ind_share) = [];
            num = num - length(ind_share);
        end
    end
end