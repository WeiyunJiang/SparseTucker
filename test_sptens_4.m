function test_sptens_4
    clear all
    clc
    
    %X = sptenrand([100 100 100 100], 50);
    %X =  creat_sptensor (0.001);
     load -ASCII nips.mat;
    %sparsity = 50/(100^4)
    tsample.nips = nips;
    i1 = nips(:,1);
    i2 = nips(:,2);
    ind_share = find(i1<=200);
    ind_share2 = find(i2(ind_share)<=200);
    nips_trunc = nips(ind_share2,:);
    
    
    ind_share3 = find(nips_trunc(:,3)<=2000);
    nips_trunc = nips_trunc(ind_share3,:);
    subs = nips_trunc(:,1:4);
    vals = nips_trunc(:,5);
    siz = [200  200  2000  17];
    tsample.X = sptensor(subs,vals,siz);
    sparsity = length(vals)/(prod(siz))
    subs = tsample.X.subs;
    vals = tsample.X.vals;
    siz = tsample.X.size;
    R = [50 50 30 13];
%     X_v = tenmat (X,[1],[2 3 4]);
%     norm_ori=norm(X_v);
       

    %%
    
    profile on

    %%
    %als
%     
     tic
     cor = tucker_als(tsample.X,R);
     toc
     err = (norm(cor)-norm(tsample.X))/norm(tsample.X)
    
%     cor_t = tensor(cor);
%     cor_v = tenmat (cor_t,[1],[2 3 4]);
%     err = norm(cor_v-X_v)/norm_ori
        %qr-p_truncated
        
        
    tic
    mine_qrp_trun = smphooi_qr_p_truncatedKron4D(tsample.X,R);
    toc
    err = (norm(mine_qrp_trun)-norm(tsample.X))/norm(tsample.X)
    
    
    
%     mine_t_qrp = tensor(mine_qrp);
%     
%     mine_v_qrp = tenmat (mine_t_qrp,[1],[2 3 4]);
%     err = norm(mine_v_qrp-X_v)/norm_ori
     %qr-p
    
    tic
    mine_qrp = smphooi_qr_p_4(tsample.X,R);
    toc
    err = (norm(mine_qrp)-norm(tsample.X))/norm(tsample.X)
%     mine_t_qrp = tensor(mine_qrp);
%     
%     mine_v_qrp = tenmat (mine_t_qrp,[1],[2 3 4]);
%     err = norm(mine_v_qrp-X_v)/norm_ori
end


 


    

    


   
   