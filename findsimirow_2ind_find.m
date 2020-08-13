function [k,j,new_sub,new_vals,ind_a,ind_b]= findsimirow_2ind_find (subs_0,vals,two_ind)
    % A is ..x2, k is the row ind
    coo_1 = subs_0(1,:);
    
    
    ind_a = coo_1(two_ind(2));
    ind_b = coo_1(two_ind(1));
    new_sub = subs_0;
    new_vals = vals;
    x2 = new_sub(:,two_ind(2));
    x3 = new_sub(:,two_ind(1));
    
    z= find(x2==ind_a);
    z_t = find(x3(z)==ind_b);
    z_b = z(z_t);
    k = new_sub(z_b,:);
    j = vals(z_b,:);
    new_sub(z_b,:) = [];
    new_vals(z_b,:) = [];
    
end