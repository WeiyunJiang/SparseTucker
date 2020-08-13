function MovieLens_trunc = truncate_MovieLens
load -ASCII ML_base1.mat

i1=ML_base1(:,1);
i2=ML_base1(:,2);
%i3=ML_base1(:,3);
%i4=ML_base1(:,4);
ind_share = find(i1<=100);
ind_share2 = find(i2(ind_share)<=200);


sparsity = length(ind_share2)/(100*200*21*24)
MovieLens_trunc = ML_base1(ind_share2,:);
end