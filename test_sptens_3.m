    clear all
    clc
    
    %X = sptenrand([100 100 100 100], 50);
    %X =  creat_sptensor_3 (0.1);
    
    %tsample = load('matmul555.mat');
    tsample = load('nell2_trun1000.mat');
    %tsample=load('lenaimagetrain.mat');
    %tsample.X = tsample.lenaimagetrain;
      subs = tsample.nell2_trun1000(:,1:3);
      vals = tsample.nell2_trun1000(:,4);
      siz = [1000 1000 1000];
     tsample.X = sptensor(subs,vals,siz);
    sparsity = length(vals)/(prod(siz))
    
 
     R = [16];
     subs = tsample.X.subs;
     vals = tsample.X.vals;
     siz=tsample.X.size;
    
    
%     X_v = tenmat (X,[1],[2 3 4]);
%     norm_ori=norm(X_v);
       

     %%

    
    profile on
     %%
    %als
    tic
    cor = tucker_als(tsample.X,R);
    toc
    err = (norm(cor)-norm(tsample.X))/norm(tsample.X)
%     cor_t = tensor(cor);
%     cor_v = tenmat (cor_t,[1],[2 3 4]);
%     err = norm(cor_v-X_v)/norm_ori

    %qr-p trunca
    tic
    mine_qrp = smphooi_qr_p_truncatedKron3D(tsample.X,R);
    toc
    err = (norm(mine_qrp)-norm(tsample.X))/norm(tsample.X)
    %parfor
    
    tic
    mine_parfor = smphooi_qr_p_parfor_truncatedKron3D(tsample.X,R);
    toc
    err = (norm(mine_parfor)-norm(tsample.X))/norm(tsample.X)
    
    %qr-p
    tic
    mine_qrp_3 = smphooi_qr_p_3(tsample.X,R);
    toc
    err = (norm(mine_qrp_3)-norm(tsample.X))/norm(tsample.X)
    
    %qr-p while
        
    tic
    mine_qrp_while = smphooi_qr_p_truncatedKron3D_while(tsample.X,R);
    toc
    err = (norm(mine_qrp_while)-norm(tsample.X))/norm(tsample.X)
    
    
%     mine_t_qrp = tensor(mine_qrp);
%     
%     mine_v_qrp = tenmat (mine_t_qrp,[1],[2 3 4]);
%     err = norm(mine_v_qrp-X_v)/norm_ori


    

    


   
   