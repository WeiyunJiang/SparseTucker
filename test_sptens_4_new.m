function test_sptens_4_new
    clear all
    clc
    
    %X = sptenrand([100 100 100 100], 50);
    %X =  creat_sptensor (0.001);
    tsample = load('Enron_trun_500_100.mat');
    %sparsity = 50/(100^4)
    
    
    subs = tsample.Enron_trun(:,1:4);
    vals = tsample.Enron_trun(:,5);
    siz = [500  500  500  100];
    tsample.X = sptensor(subs,vals,siz);
    sparsity = length(vals)/(prod(siz))
    subs = tsample.X.subs;
    vals = tsample.X.vals;
    siz = tsample.X.size;
    R = 20
%     X_v = tenmat (X,[1],[2 3 4]);
%     norm_ori=norm(X_v);
       

    %%
    
    profile on

    %%
    %als
%     
%      tic
%      cor = tucker_als(tsample.X,R);
%      toc
%      err = (norm(cor)-norm(tsample.X))/norm(tsample.X)
%      [sparsity_core1,sparsity_fac1]=sparsity_core_factor(cor,4)
%     cor_t = tensor(cor);
%     cor_v = tenmat (cor_t,[1],[2 3 4]);
%     err = norm(cor_v-X_v)/norm_ori
        %qr-p_truncated
        
        
    tic
    mine_qrp_trun = smphooi_qr_p_truncatedKron4D(tsample.X,R);
    toc
    err = (norm(mine_qrp_trun)-norm(tsample.X))/norm(tsample.X)
    [sparsity_core2,sparsity_fac2]=sparsity_core_factor(mine_qrp_trun,4)
    
    
%     mine_t_qrp = tensor(mine_qrp);
%     
%     mine_v_qrp = tenmat (mine_t_qrp,[1],[2 3 4]);
%     err = norm(mine_v_qrp-X_v)/norm_ori
     %qr-p
    
    tic
    mine_qrp = smphooi_qr_p_4(tsample.X,R);
    toc
    err = (norm(mine_qrp)-norm(tsample.X))/norm(tsample.X)
    [sparsity_core3,sparsity_fac3]=sparsity_core_factor(mine_qrp,4)
%     mine_t_qrp = tensor(mine_qrp);
%     
%     mine_v_qrp = tenmat (mine_t_qrp,[1],[2 3 4]);
%     err = norm(mine_v_qrp-X_v)/norm_ori

end


 


    

    


   
   