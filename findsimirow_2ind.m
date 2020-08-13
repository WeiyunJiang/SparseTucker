function [k,j,new_sub,new_vals,ind_a,ind_b]= findsimirow_2ind (subs_0,vals,two_ind)
    % A is ..x2, k is the row ind
    coo_1 = subs_0(1,:);
    height = size(subs_0,1);
    
    ind_a = coo_1(two_ind(2));
    ind_b = coo_1(two_ind(1));
    new_sub = subs_0;
    new_vals = vals;
    k_ind = 1;
    
    k = zeros(1,3);
    j=0;
    del_ind =0;
    for i = 2:height
        new_line =subs_0(i,:);
        if ind_a == new_line(two_ind(2)) && ind_b == new_line(two_ind(1))
            k(k_ind,:) = new_line;
            j(k_ind) = vals(i);
            k_ind = k_ind+1;
            new_sub(i+del_ind,:) =[]; 
            new_vals(i+del_ind) = [];
            del_ind = del_ind-1;
        end
    end
    new_sub(1,:) = [];
    new_vals(1) = [];
    k(k_ind,:) = coo_1;
    j(k_ind) = vals(1);
    j=j';
end