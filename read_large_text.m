function  read_large_text

clear 
clc
clear all

load -ASCII nell-1.mat
subs = chicago_crime_comm(:,1:4);
vals = chicago_crime_comm(:,5);
siz = [6186 24  77  32];
X = sptensor(subs,vals,siz);

tic 
X_my = tensor(smphooi_qr_p_truncatedKron4D(X,16));
toc

err = (norm(X_my)-norm(X))/norm(X)

tic
als = tensor(tucker_als(X,16));
toc 
err = (norm(als)-norm(X))/norm(X)


end