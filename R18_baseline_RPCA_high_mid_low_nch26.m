% same as R12, but only fewer channels, excluding 26,27,29,30.
% (1) load graphs, edge = lowVOT ...looking at increases and decreases
% (2) load graphs, edge = midVOT ...looking at increases and decreases
% (3) load graphs, edge = highVOT ...looking at increases and decreases

chlist = [1:25 28];
%% (1) load graphs, edge = lowVOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Robust PCA in time%%%%%%%%%%%%%%%%%%

miVOT123 = load('mi_matrix_VOT123_K4.mat'); %with mem 2
miVOT456 = load('mi_matrix_VOT456_K4.mat'); 
np = 22; %number of people 21 undergrads, 1 aphasia
nwin = 43;
nch = length(chlist);
ntim = [1 51 101:5:301];
% creating the adjacency matrices, with edge = VOT456- VOT123, only
% positive increases diff_mid_lowVOT
diff_lowVOT = zeros(np,nch,nch,nwin-2);
diff_mid_lowVOT = zeros(np,nch,nch,nwin-2);

for ipat = 1:np
    mi123 = miVOT123.mi_matrix_VOT123_K4{ipat,1}; %contains 43 time windows
    mi456 = miVOT456.mi_matrix_VOT456_K4{ipat,1};
    for iwin = 3:nwin
       mat1 = squeeze(mi456{iwin,1} - mi456{1,1}); 
       mat2 = squeeze(mi123{iwin,1} - mi123{1,1}); 
       diff_lowVOT(ipat,:,:,iwin-2)= squeeze(mat2(chlist,chlist));
       mat =  squeeze(mat1-mat2);
       diff_mid_lowVOT(ipat,:,:,iwin-2) =mat(chlist, chlist);
    end

end

nwin = nwin-2;
energy_Stime = zeros(np,nwin);
energy_Xtime = zeros(np,nwin);
S_pat_ml_lm = cell(1,np);
L_pat_ml_lm = cell(1,np);
parfor ipat = 1:np
    tic
    X = squeeze(diff_mid_lowVOT(ipat,:,:,:));
    [Lx,Sx] = RobustPCA_time(X);
    Y = squeeze(diff_lowVOT(ipat,:,:,:));
    [Ly,Sy] = RobustPCA_time(Y);
    S_pat_ml_lm{1,ipat}=Sx;
    L_pat_ml_lm{1,ipat}=Lx;
    for iwin = 1:nwin
        [U1,S1,V1] = eig(squeeze(Sx(:,:,iwin)));
        energy_S_mid_low(ipat,iwin) = sum(abs(diag(S1)));
        [U2,S2,V2] = eig(squeeze(Sy(:,:,iwin)));
        energy_S_low(ipat,iwin) = sum(abs(diag(S2)));
    end
    toc
end


X1 = energy_S_mid_low(:,1:end);
X2 = energy_S_low(:,1:end);
figure; 
 plot(X1(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X1(end,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S mid-low');

figure; 
 plot(X2(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X2(end,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S just low');

[coeff, score, latent, tsquared, explained, mu] = pca(X1(:,1:end));%20
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('baseline norm Edge = mid-low');

[coeff, score, latent, tsquared, explained, mu] = pca(X2(:,1:end));%20
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('baseline norm Edge = low');

% reconstruct with 1st PC
[coeff, score, latent, tsquared, explained, mu] = pca(X1(:,1:end));%20
Xhat = score(:,1:3)*coeff(:,1:3)' + mu;
 figure;
 plot(Xhat(1:end-1,:)','-b');hold on;
plot(Xhat(end,:)','-r');hold on;

%% (2) load graphs, edge = midVOT - highVOT ...looking at increases and decreases
%%%% load new graphs, edge = midVOT - highVOT ...looking at increases only
miVOT789 = load('mi_matrix_VOT789_K4.mat'); %with mem 2
miVOT456 = load('mi_matrix_VOT456_K4.mat'); 
np = 22; %number of people 21 undergrads, 1 aphasia
nwin = 43;
nch = 26;
ntim = [1 51 101:5:301];
% creating the adjacency matrices, with edge = VOT456- VOT789, only
% positive increases diff_mid_highVOT
diff_midVOT = zeros(np,nch,nch,nwin-2);
diff_mid_highVOT = zeros(np,nch,nch,nwin-2);

for ipat = 1:np
    mi789 = miVOT789.mi_matrix_VOT789_K4{ipat,1}; %contains 43 time windows
    mi456 = miVOT456.mi_matrix_VOT456_K4{ipat,1};
    for iwin = 3:nwin
       mat1 = squeeze(mi456{iwin,1} -mi456{1,1} ); 
       mat2 = squeeze(mi789{iwin,1}- mi789{1,1}); 
       diff_midVOT(ipat,:,:,iwin-2)= squeeze(mat1(chlist,chlist));
       mat = squeeze(mat1-mat2);
       diff_mid_highVOT(ipat,:,:,iwin-2) = mat(chlist,chlist);
    end

end

nwin = nwin-2;
energy_S_mid_high = zeros(np,nwin);
energy_S_mid = zeros(np,nwin);
S_pat_mh_hm = cell(1,np);
L_pat_mh_hm = cell(1,np);
parfor ipat = 1:np
    tic
    X = squeeze(diff_mid_highVOT(ipat,:,:,:));
    [Lx,Sx] = RobustPCA_time(X);
    Y = squeeze(diff_midVOT(ipat,:,:,:));
    [Ly,Sy] = RobustPCA_time(Y);
    S_pat_mh_hm{1,ipat}=Sx;
    L_pat_mh_hm{1,ipat}=Lx;
    for iwin = 1:nwin
        [U1,S1,V1] = eig(squeeze(Sx(:,:,iwin)));
        energy_S_mid_high(ipat,iwin) = sum(abs(diag(S1)));
        [U2,S2,V2] = eig(squeeze(Sy(:,:,iwin)));
        energy_S_mid(ipat,iwin) = sum(abs(diag(S2)));
    end
    toc
end


X1 = energy_S_mid_high(:,1:end);
X2 = energy_S_mid(:,1:end);
figure; 
 plot(X1(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X1(end,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S mid-high');

figure; 
 plot(X2(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X2(end,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S just mid');

[coeff, score, latent, tsquared, explained, mu] = pca(X1(:,1:end));%20
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('baseline norm Edge = mid-high');

[coeff, score, latent, tsquared, explained, mu] = pca(X2(:,1:end));%20
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('baseline norm Edge = mid');

% reconstruct with 1st PC
[coeff, score, latent, tsquared, explained, mu] = pca(X1(:,1:end));%20
Xhat = score(:,1:3)*coeff(:,1:3)' + mu;
 figure;
 plot(Xhat(1:end-1,:)','-b');hold on;
plot(Xhat(end,:)','-r');hold on;

%% (3) load graphs, edge = lowVOT - highVOT ...looking at increases and decreases
%%%% load new graphs, edge = lowVOT - highVOT ...looking at increases only
miVOT123 = load('mi_matrix_VOT123_K4.mat'); %with mem 2
miVOT789 = load('mi_matrix_VOT789_K4.mat'); 
np = 22; %number of people end undergrads, 1 aphasia
nwin = 43;
nch = 26;
ntim = [1 51 101:5:301];
% creating the adjacency matrices, with edge = VOT456- VOT789, only
% positive increases diff_low_highVOT
diff_highVOT = zeros(np,nch,nch,nwin-2);
diff_low_highVOT = zeros(np,nch,nch,nwin-2);

for ipat = 1:np
    mi789 = miVOT789.mi_matrix_VOT789_K4{ipat,1}; %contains 43 time windows
    mi123 = miVOT123.mi_matrix_VOT123_K4{ipat,1};
    for iwin = 3:nwin
       mat1 = squeeze(mi123{iwin,1}- mi123{1,1}); 
       mat2 = squeeze(mi789{iwin,1}- mi789{1,1}); 
       diff_highVOT(ipat,:,:,iwin-2)= squeeze(mat2(chlist,chlist));
       mat = squeeze(mat1-mat2);
       diff_low_highVOT(ipat,:,:,iwin-2) = mat(chlist,chlist);
    end

end

nwin = nwin-2;
energy_S_low_high = zeros(np,nwin);
energy_S_high = zeros(np,nwin);
S_pat_lh_hl = cell(1,np);
L_pat_lh_hl = cell(1,np);
parfor ipat = 1:np
    tic
    X = squeeze(diff_low_highVOT(ipat,:,:,:));
    [Lx,Sx] = RobustPCA_time(X);
    Y = squeeze(diff_highVOT(ipat,:,:,:));
    [Ly,Sy] = RobustPCA_time(Y);
    S_pat_lh_hl{1,ipat}=Sx;
    L_pat_lh_hl{1,ipat}=Lx;
    for iwin = 1:nwin
        [U1,S1,V1] = eig(squeeze(Sx(:,:,iwin)));
        energy_S_low_high(ipat,iwin) = sum(abs(diag(S1)));
        [U2,S2,V2] = eig(squeeze(Sy(:,:,iwin)));
        energy_S_high(ipat,iwin) = sum(abs(diag(S2)));
    end
    toc
end


X1 = energy_S_low_high(:,1:end);
X2 = energy_S_high(:,1:end);
figure; 
 plot(X1(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X1(end,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S low-high');

figure; 
 plot(X2(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X2(end,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S just high');

[coeff, score, latent, tsquared, explained, mu] = pca(X1(:,1:end));%20
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('baseline norm Edge = low-high');

[coeff, score, latent, tsquared, explained, mu] = pca(X2(:,1:end));%20
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('baseline norm Edge = high');

% reconstruct with 1st PC
[coeff, score, latent, tsquared, explained, mu] = pca(X1(:,1:end));%20
Xhat = score(:,1:3)*coeff(:,1:3)' + mu;
 figure;
 plot(Xhat(1:end-1,:)','-b');hold on;
plot(Xhat(end,:)','-r');hold on;

%% save the L and S matrices for all 3 conditions.
BS_RpcaN26.S_pat_ml_lm = S_pat_ml_lm;
BS_RpcaN26.L_pat_ml_lm = L_pat_ml_lm;
BS_RpcaN26.S_pat_mh_hm = S_pat_mh_hm;
BS_RpcaN26.L_pat_mh_hm = L_pat_mh_hm;
BS_RpcaN26.S_pat_lh_hl = S_pat_lh_hl;
BS_RpcaN26.L_pat_lh_hl = L_pat_lh_hl;
    
save('BS_RpcaN26','BS_RpcaN26','-v7.3');
%% this program compares baseline normalized matrices using edges mid-low, etc.
load('BS_RpcaN26');
energy_S_low_high = zeros(np,nwin);
energy_L_low_high = zeros(np,1);
for ipat = 1:np
    Sx = BS_RpcaN26.S_pat_lh_hl{1,ipat};
    Lx = BS_RpcaN26.L_pat_lh_hl{1,ipat};
    [U2,S2,V2] = eig(squeeze(Lx(:,:)));
    energy_L_low_high(ipat,1) = sum(abs(diag(S2)));
    
    for iwin = 1:nwin
        [U1,S1,V1] = eig(squeeze(Sx(:,:,iwin)));
        energy_S_low_high(ipat,iwin) = sum(abs(diag(S1)));
    end  
end

X1 = energy_S_low_high(:,1:end);
ntim = [101:5:301];

figure; 
plot(X1(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X1(end,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S low-high');

figure; plot(ntim,mean(X1(1:end,:)),'*-k','LineWidth',1);hold on;
plot(ntim, X1(end,:),'-r*','LineWidth',1);
[ph,msg]=jbfill(ntim,(mean(X1(1:end,:))+std(X1(1:end,:))),(mean(X1(1:end,:))-std(X1(1:end,:))),'b','b',1,0.3);

stats = zeros(1,nwin);
stats_lh = zeros(1,np);
for ipat = 1:np
    pat_list = 1:np;
    pat_list(ipat)=[];
    for iwin = 1:nwin
        [h,p, t, df] =  ttestch(mean(X1(pat_list,iwin)), std(X1(pat_list,iwin)), X1(ipat,iwin), np-1, 0.05);
        stats(iwin) = p;
        
    end
    stats = fdr(stats);
    fprintf('low-high %d, numb of p<0.05 = %d\n',ipat, sum(stats<0.05));
    stats_lh(ipat) = sum(stats<0.05);
end

%%%%%%%%%%%%%%
stats_ml = zeros(1,np);
energy_S_mid_low = zeros(np,nwin);
energy_L_mid_low = zeros(np,1);
for ipat = 1:np
    Sx = BS_RpcaN26.S_pat_ml_lm{1,ipat};
    Lx = BS_RpcaN26.L_pat_ml_lm{1,ipat};
    [U2,S2,V2] = eig(squeeze(Lx(:,:)));
    energy_L_mid_low(ipat,1) = sum(abs(diag(S2)));
    
    for iwin = 1:nwin
        [U1,S1,V1] = eig(squeeze(Sx(:,:,iwin)));
        energy_S_mid_low(ipat,iwin) = sum(abs(diag(S1)));
    end  
end

X1 = energy_S_mid_low(:,1:end);
ntim = [101:5:301];

figure; 
plot(X1(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X1(end,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S mid_low');

figure; plot(ntim,mean(X1(1:end,:)),'*-k','LineWidth',1);hold on;
plot(ntim, X1(end,:),'-r*','LineWidth',1);
[ph,msg]=jbfill(ntim,(mean(X1(1:end,:))+std(X1(1:end,:))),(mean(X1(1:end,:))-std(X1(1:end,:))),'b','b',1,0.3);

stats = zeros(1,nwin);

for ipat = 1:np
    pat_list = 1:np;
    pat_list(ipat)=[];
    for iwin = 1:nwin
        [h,p, t, df] =  ttestch(mean(X1(pat_list,iwin)), std(X1(pat_list,iwin)), X1(ipat,iwin), np-1, 0.05);
        stats(iwin) = p;
        
    end
    stats = fdr(stats);
    fprintf('mid_low %d, numb of p<0.05 = %d\n',ipat, sum(stats<0.05));
    stats_ml(ipat) = sum(stats<0.05);
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stats_mh = zeros(1,np);
energy_S_mid_high = zeros(np,nwin);
energy_L_mid_high = zeros(np,1);
for ipat = 1:np
    Sx = BS_RpcaN26.S_pat_mh_hm{1,ipat};
    Lx = BS_RpcaN26.L_pat_mh_hm{1,ipat};
    [U2,S2,V2] = eig(squeeze(Lx(:,:)));
    energy_L_mid_high(ipat,1) = sum(abs(diag(S2)));
    
    for iwin = 1:nwin
        [U1,S1,V1] = eig(squeeze(Sx(:,:,iwin)));
        energy_S_mid_high(ipat,iwin) = sum(abs(diag(S1)));
    end  
end

X1 = energy_S_mid_high(:,1:end);
ntim = [101:5:301];

figure; 
plot(X1(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X1(end,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S mid_high');

figure; plot(ntim,mean(X1(1:end,:)),'*-k','LineWidth',1);hold on;
plot(ntim, X1(end,:),'-r*','LineWidth',1);
[ph,msg]=jbfill(ntim,(mean(X1(1:end,:))+std(X1(1:end,:))),(mean(X1(1:end,:))-std(X1(1:end,:))),'b','b',1,0.3);

stats = zeros(1,nwin);

for ipat = 1:np
    pat_list = 1:np;
    pat_list(ipat)=[];
    for iwin = 1:nwin
        [h,p, t, df] =  ttestch(mean(X1(pat_list,iwin)), std(X1(pat_list,iwin)), X1(ipat,iwin), np-1, 0.05);
        stats(iwin) = p;
        
    end
    stats = fdr(stats);
    fprintf('mid_high %d, numb of p<0.05 = %d\n',ipat, sum(stats<0.05));
    stats_mh(ipat) = sum(stats<0.05);
end
    

%% plotting figures with statistics
figure;
subplot(3,1,1);hold on;
bar(stats_lh);
plot(1:np,stats_lh,'ko')
subplot(3,1,2);hold on;
bar(stats_ml);
plot(1:np,stats_ml,'ko')
subplot(3,1,3);hold on;
bar(stats_mh);
plot(1:np,stats_mh,'ko')
