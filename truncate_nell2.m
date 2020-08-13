function nell2_trun = truncate_nell2
sample =load('nell2.mat');
tru_n = 2000;
i1=sample.nell2(:,1);
i2=sample.nell2(:,2);
i3=sample.nell2(:,3);
ind_share = find(i1<=tru_n);
ind_share2 = find(i2(ind_share)<=tru_n);
nell2_trun = sample.nell2(ind_share2,:);
ind_share3 = find(nell2_trun(:,3)<=tru_n);
sparsity = length(ind_share3)/tru_n^3
nell2_trun = nell2_trun(ind_share3,:);
end
