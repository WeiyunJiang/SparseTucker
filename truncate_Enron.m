function Enron_trun = truncate_Enron
load -ASCII enron.mat;
i1=enron(:,1);
i2=enron(:,2);
i3=enron(:,3);
i4=enron(:,4);
tru_n=1000;
ind_share = find(i1<=tru_n);
ind_share2 = find(i2(ind_share)<=tru_n);
Enron_trun = enron(ind_share2,:);
ind_share3 = find(Enron_trun(:,3)<=tru_n);
Enron_trun = Enron_trun(ind_share3,:);
ind_share4 = find(Enron_trun(:,4)<=100);
Enron_trun = Enron_trun(ind_share4,:);
sparsity = length(ind_share4)/tru_n^4

end