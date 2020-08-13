function nell1_trun = truncate_nell1
sample =load('nell1.mat');
tru_n = 2000;
i1=sample.nell1(:,1);
i2=sample.nell1(:,2);
i3=sample.nell1(:,3);
ind_share = find(i1<=tru_n);
ind_share2 = find(i2(ind_share)<=tru_n);
nell1_trun = sample.nell2(ind_share2,:);
ind_share3 = find(nell1_trun(:,3)<=tru_n);
sparsity = length(ind_share3)/tru_n^3
nell1_trun = nell1_trun(ind_share3,:);
end