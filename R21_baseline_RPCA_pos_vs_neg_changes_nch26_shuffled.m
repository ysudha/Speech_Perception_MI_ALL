% same as R19 fewer channels. - but shuffling the matrices.
%do incerases in MI in (mid-high) or (mid-low) or (high-low) matter more or do
% decreases matter more?
% (1) load graphs, edge = lowVOT ...looking at increases, decreases
% (2) load graphs, edge = midVOT ...looking at increases, decreases
% (3) load graphs, edge = highVOT ...looking at increases, decreases

chlist = [1:25 28];%[1 3 4 8 12 13 2 6 7 11 15 16];%[1:25 28];
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
       mat = squeeze(mat1-mat2);
       diff_mid_lowVOT(ipat,:,:,iwin-2) = mat(chlist,chlist);
    end

end

nwin = nwin-2;
energy_S_mid_gt_low = zeros(np,nwin);
energy_S_mid_less_low = zeros(np,nwin);
Q_S_mid_gt_low = zeros(np,nwin);
Q_S_mid_less_low = zeros(np,nwin);
parfor ipat = 1:np
    tic
    X = squeeze(diff_mid_lowVOT(ipat,:,:,:));
    X(X<0)=0;
    [Lx,Sx] = RobustPCA_time(X);
    Y = squeeze(diff_mid_lowVOT(ipat,:,:,:));
    Y(Y>0)=0;
    Y=Y.*(-1);
    [Ly,Sy] = RobustPCA_time(Y);
    for iwin = 1:nwin
        [U1,S1,V1] = eig(squeeze(Sx(:,:,iwin)));
        energy_S_mid_gt_low(ipat,iwin) = sum(abs(diag(S1)));
        [partition1,Q_S_mid_gt_low(ipat,iwin)]=community_louvain(squeeze(X(:,:,iwin)),1,randperm(nch));
        [U2,S2,V2] = eig(squeeze(Sy(:,:,iwin)));
        energy_S_mid_less_low(ipat,iwin) = sum(abs(diag(S2)));
        [partition1,Q_S_mid_less_low(ipat,iwin)]=community_louvain(squeeze(Y(:,:,iwin)),1,randperm(nch));

    end
    toc
end


X1 = energy_S_mid_gt_low(:,1:end);
X2 = energy_S_mid_less_low(:,1:end);
figure; 
 plot(X1(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X1(end,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S mid-low >0');

figure; 
 plot(X2(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X2(end,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S mid-low<0');

[coeff, score, latent, tsquared, explained, mu] = pca(X1(:,1:end));%20
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('baseline norm Edge = mid-low >0');

[coeff, score, latent, tsquared, explained, mu] = pca(X2(:,1:end));%20
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('baseline norm Edge = mid-low<0');


X1 = Q_S_mid_gt_low(:,1:end);
X2 = Q_S_mid_less_low(:,1:end);
figure; 
 plot(X1(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X1(end,:)','-r*','LineWidth',4);hold on;
title('Q of mid-low >0');

figure; 
 plot(X2(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X2(end,:)','-r*','LineWidth',4);hold on;
title('Q of mid-low<0');

[coeff, score, latent, tsquared, explained, mu] = pca(X1(:,1:end));%20
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('Q = mid-low >0');

[coeff, score, latent, tsquared, explained, mu] = pca(X2(:,1:end));%20
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('Q = mid-low<0');

% similarity between Qvariations and energy of S
figure;hold on;
for i = np
    scatter(energy_S_mid_gt_low(i,:),Q_S_mid_gt_low(i,:));
    %corr(energy_S_mid_gt_low(i,:)',Q_S_mid_gt_low(i,:)')
end
% % patient similarity matrix
% pat_sim_mid_gt_low = corr(X1',X1','Type','Spearman');
% pat_sim_mid_gt_low = pat_sim_mid_gt_low - diag(diag(pat_sim_mid_gt_low));
% figure; imagesc(pat_sim_mid_gt_low); colormap jet;
% 
% pat_sim_mid_lt_low = corr(X2',X2','Type','Spearman');
% pat_sim_mid_lt_low = pat_sim_mid_lt_low - diag(diag(pat_sim_mid_lt_low));
% figure; imagesc(pat_sim_mid_lt_low); colormap jet;
%% (2) load graphs, edge = midVOT - highVOT ...looking at increases and decreases
%%%% load new graphs, edge = midVOT - highVOT ...looking at increases only
miVOT789 = load('mi_matrix_VOT789_K4.mat'); %with mem 2
miVOT456 = load('mi_matrix_VOT456_K4.mat'); 
np = 22; %number of people 21 undergrads, 1 aphasia
nwin = 43;
nch = length(chlist);
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
       mat = squeeze(mat1-mat2);
       diff_mid_highVOT(ipat,:,:,iwin-2) = mat(chlist,chlist);
    end

end

nwin = nwin-2;
energy_S_mid_gt_high = zeros(np,nwin);
energy_S_mid_less_high = zeros(np,nwin);
parfor ipat = 1:np
    tic
    X = squeeze(diff_mid_highVOT(ipat,:,:,:));
    X(X<0)=0;
    [Lx,Sx] = RobustPCA_time(X);
    Y = squeeze(diff_mid_highVOT(ipat,:,:,:));
    Y(Y>0)=0;
    [Ly,Sy] = RobustPCA_time(Y);
    for iwin = 1:nwin
        [U1,S1,V1] = eig(squeeze(Sx(:,:,iwin)));
        energy_S_mid_gt_high(ipat,iwin) = sum(abs(diag(S1)));
        [U2,S2,V2] = eig(squeeze(Sy(:,:,iwin)));
        energy_S_mid_less_high(ipat,iwin) = sum(abs(diag(S2)));
    end
    toc
end


X1 = energy_S_mid_gt_high(:,1:end);
X2 = energy_S_mid_less_high(:,1:end);
figure; 
 plot(X1(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X1(end,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S mid-high pos');

figure; 
 plot(X2(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X2(end,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S mid-high neg');

[coeff, score, latent, tsquared, explained, mu] = pca(X1(:,1:end));%20
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('baseline norm Edge = mid-high pos');

[coeff, score, latent, tsquared, explained, mu] = pca(X2(:,1:end));%20
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('baseline norm Edge = mid-high neg');

% patient similarity matrix
% pat_sim_mid_gt_high = corr(X1',X1','Type','Spearman');
% pat_sim_mid_gt_high = pat_sim_mid_gt_high - diag(diag(pat_sim_mid_gt_high));
% pat_sim_mid_gt_high(pat_sim_mid_gt_high<0.3)=0;
% figure; imagesc(pat_sim_mid_gt_high); colormap jet;
% 
% pat_sim_mid_lt_high = corr(X2',X2','Type','Spearman');
% pat_sim_mid_lt_high = pat_sim_mid_lt_high - diag(diag(pat_sim_mid_lt_high));
% pat_sim_mid_lt_high(pat_sim_mid_lt_high<0.3)=0;
% figure; imagesc(pat_sim_mid_lt_high); colormap jet;
%% (3) load graphs, edge = lowVOT - highVOT ...looking at increases and decreases
%%%% load new graphs, edge = lowVOT - highVOT ...looking at increases only
miVOT123 = load('mi_matrix_VOT123_K4.mat'); %with mem 2
miVOT789 = load('mi_matrix_VOT789_K4.mat'); 
np = 22; %number of people end undergrads, 1 aphasia
nwin = 43;
nch = length(chlist);
ntim = [1 51 101:5:301];
% creating the adjacency matrices, with edge = VOT456- VOT789, only
% positive increases diff_low_highVOT
diff_low_highVOT = zeros(np,nch,nch,nwin-2);

for ipat = 1:np
    mi789 = miVOT789.mi_matrix_VOT789_K4{ipat,1}; %contains 43 time windows
    mi123 = miVOT123.mi_matrix_VOT123_K4{ipat,1};
    for iwin = 3:nwin
       mat1 = squeeze(mi123{iwin,1}- mi123{1,1}); 
       mat2 = squeeze(mi789{iwin,1}- mi789{1,1}); 
       mat = squeeze(mat1-mat2);
       diff_low_highVOT(ipat,:,:,iwin-2) = mat(chlist,chlist);
    end

end

nwin = nwin-2;
energy_S_low_gt_high = zeros(np,nwin);
energy_S_low_less_high = zeros(np,nwin);
S_pat = cell(1,np);
L_pat = cell(1,np);
parfor ipat = 1:np
    tic
    X = squeeze(diff_low_highVOT(ipat,:,:,:));
    [Lx,Sx] = RobustPCA_time(X);
    X(X<0)=0;
    Y = squeeze(diff_low_highVOT(ipat,:,:,:));
    Y(Y>0)=0;
    [Ly,Sy] = RobustPCA_time(Y);
    for iwin = 1:nwin
        [U1,S1,V1] = eig(squeeze(Sx(:,:,iwin)));
        energy_S_low_gt_high(ipat,iwin) = sum(abs(diag(S1)));
        [U2,S2,V2] = eig(squeeze(Sy(:,:,iwin)));
        energy_S_low_less_high(ipat,iwin) = sum(abs(diag(S2)));
    end
    toc
end


X1 = energy_S_low_gt_high(:,1:end);
X2 = energy_S_low_less_high(:,1:end);
figure; 
 plot(X1(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X1(end,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S low-high>0');

figure; 
 plot(X2(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X2(end,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S low-high<0');

[coeff, score, latent, tsquared, explained, mu] = pca(X1(:,1:end));%20
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('baseline norm Edge = low-high>0');

[coeff, score, latent, tsquared, explained, mu] = pca(X2(:,1:end));%20
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('baseline norm Edge = low-high<0');

%%%%%%%%%%% saving the data%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RPCA_pos_negN26.energy_S_mid_gt_high = energy_S_mid_gt_high;
RPCA_pos_negN26.energy_S_mid_less_high = energy_S_mid_less_high;
RPCA_pos_negN26.energy_S_mid_gt_low=energy_S_mid_gt_low;
RPCA_pos_negN26.energy_S_mid_less_low = energy_S_mid_less_low;
RPCA_pos_negN26.energy_S_low_gt_high = energy_S_low_gt_high;
RPCA_pos_negN26.energy_S_low_less_high = energy_S_low_less_high;

save('RPCA_pos_negN26','RPCA_pos_negN26','-v7.3');
%% plotting code
load('RPCA_pos_negN26');

%% plotting code (1) 
X1 = RPCA_pos_negN26.energy_S_mid_gt_high(:,1:end);
X2 = RPCA_pos_negN26.energy_S_mid_less_high(:,1:end);
ntim = ([101:5:301]-101)*2;

figure; 
plot(X1(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X1(end,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S mid_gt_high');

[coeff, score, latent, tsquared, explained, mu] = pca(X1(:,1:end));%20
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('baseline norm Edge = mid_gt_high');

figure; plot(ntim,mean(X1(1:end,:)),'*-k','LineWidth',1);hold on;
plot(ntim, X1(end,:),'-r*','LineWidth',2);
[ph,msg]=jbfill(ntim,(mean(X1(1:end,:))+std(X1(1:end,:))),(mean(X1(1:end,:))-std(X1(1:end,:))),'b','b',1,0.3);
title('baseline norm Energy of S mid_gt_high');
xlabel('Time(ms)');

stats = zeros(1,nwin);

for ipat = 1:np
    pat_list = 1:np;
    pat_list(ipat)=[];
    for iwin = 1:nwin
        [h,p, t, df] =  ttestch(mean(X1(pat_list,iwin)), std(X1(pat_list,iwin)), X1(ipat,iwin), np-1, 0.05);
        stats(iwin) = p;
        
    end
    stats = fdr(stats);
    fprintf('mid_gt_high %d, numb of p<0.05 = %d\n',ipat, sum(stats<0.05));

end
%%%%%%%%%
figure; 
plot(X2(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X2(end,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S mid_less_high');

[coeff, score, latent, tsquared, explained, mu] = pca(X2(:,1:end));%20
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('baseline norm Edge = mid_less_high');

figure; plot(ntim,mean(X2(1:end,:)),'*-k','LineWidth',1);hold on;
plot(ntim, X2(end,:),'-r*','LineWidth',2);
[ph,msg]=jbfill(ntim,(mean(X2(1:end,:))+std(X2(1:end,:))),(mean(X2(1:end,:))-std(X2(1:end,:))),'b','b',1,0.3);
title('baseline norm Energy of S mid_less_high');
xlabel('Time(ms)');

stats = zeros(1,nwin);

for ipat = 1:np
    pat_list = 1:np;
    pat_list(ipat)=[];
    for iwin = 1:nwin
        [h,p, t, df] =  ttestch(mean(X2(pat_list,iwin)), std(X2(pat_list,iwin)), X2(ipat,iwin), np-1, 0.05);
        stats(iwin) = p;
        
    end
    stats = fdr(stats);
    fprintf('mid_less_high %d, numb of p<0.05 = %d\n',ipat, sum(stats<0.05));

end
%%%%%%%%%%%%%%
%% plotting (2)

X1 = RPCA_pos_negN26.energy_S_mid_less_low(:,1:end);
X2 = RPCA_pos_negN26.energy_S_mid_gt_low(:,1:end);
ntim = ([101:5:301]-101)*2;

figure; 
plot(X1(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X1(end,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S mid_less_low');

[coeff, score, latent, tsquared, explained, mu] = pca(X1(:,1:end));%20
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('baseline norm Edge = mid_less_low');

figure; plot(ntim,mean(X1(1:end,:)),'*-k','LineWidth',1);hold on;
plot(ntim, X1(end,:),'-r*','LineWidth',2);
[ph,msg]=jbfill(ntim,(mean(X1(1:end,:))+std(X1(1:end,:))),(mean(X1(1:end,:))-std(X1(1:end,:))),'b','b',1,0.3);
title('baseline norm Energy of S mid_less_low');
xlabel('Time(ms)');

stats = zeros(1,nwin);

for ipat = 1:np
    pat_list = 1:np;
    pat_list(ipat)=[];
    for iwin = 1:nwin
        [h,p, t, df] =  ttestch(mean(X1(pat_list,iwin)), std(X1(pat_list,iwin)), X1(ipat,iwin), np-1, 0.05);
        stats(iwin) = p;
        
    end
    stats = fdr(stats);
    fprintf('mid_less_low %d, numb of p<0.05 = %d\n',ipat, sum(stats<0.05));

end
%%%%%%%%%
figure; 
plot(X2(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X2(end,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S mid_gt_low');

[coeff, score, latent, tsquared, explained, mu] = pca(X2(:,1:end));%20
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('baseline norm Edge = mid_gt_low');

figure; plot(ntim,mean(X2(1:end,:)),'*-k','LineWidth',1);hold on;
plot(ntim, X2(end,:),'-r*','LineWidth',2);
[ph,msg]=jbfill(ntim,(mean(X2(1:end,:))+std(X2(1:end,:))),(mean(X2(1:end,:))-std(X2(1:end,:))),'b','b',1,0.3);
title('baseline norm Energy of S mid_gt_low');
xlabel('Time(ms)');

stats = zeros(1,nwin);

for ipat = 1:np
    pat_list = 1:np;
    pat_list(ipat)=[];
    for iwin = 1:nwin
        [h,p, t, df] =  ttestch(mean(X2(pat_list,iwin)), std(X2(pat_list,iwin)), X2(ipat,iwin), np-1, 0.05);
        stats(iwin) = p;
        
    end
    stats = fdr(stats);
    fprintf('mid_gt_low %d, numb of p<0.05 = %d\n',ipat, sum(stats<0.05));

end

%*****************
%% plotting (3)

X1 = RPCA_pos_negN26.energy_S_low_less_high(:,1:end);
X2 = RPCA_pos_negN26.energy_S_low_gt_high(:,1:end);
ntim = ([101:5:301]-101)*2;

figure; 
plot(X1(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X1(end,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S low_less_high');

[coeff, score, latent, tsquared, explained, mu] = pca(X1(:,1:end));%20
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('baseline norm Edge = low_less_high');

figure; plot(ntim,mean(X1(1:end,:)),'*-k','LineWidth',1);hold on;
plot(ntim, X1(end,:),'-r*','LineWidth',2);
[ph,msg]=jbfill(ntim,(mean(X1(1:end,:))+std(X1(1:end,:))),(mean(X1(1:end,:))-std(X1(1:end,:))),'b','b',1,0.3);
title('baseline norm Energy of S low_less_high');
xlabel('Time(ms)');

stats = zeros(1,nwin);

for ipat = 1:np
    pat_list = 1:np;
    pat_list(ipat)=[];
    for iwin = 1:nwin
        [h,p, t, df] =  ttestch(mean(X1(pat_list,iwin)), std(X1(pat_list,iwin)), X1(ipat,iwin), np-1, 0.05);
        stats(iwin) = p;
        
    end
    stats = fdr(stats);
    fprintf('low_less_high %d, numb of p<0.05 = %d\n',ipat, sum(stats<0.05));

end
%%%%%%%%%
figure; 
plot(X2(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X2(end,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S low_gt_high');

[coeff, score, latent, tsquared, explained, mu] = pca(X2(:,1:end));%20
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('baseline norm Edge = low_gt_high');

figure; plot(ntim,mean(X2(1:end,:)),'*-k','LineWidth',1);hold on;
plot(ntim, X2(end,:),'-r*','LineWidth',2);
[ph,msg]=jbfill(ntim,(mean(X2(1:end,:))+std(X2(1:end,:))),(mean(X2(1:end,:))-std(X2(1:end,:))),'b','b',1,0.3);
title('baseline norm Energy of S low_gt_high');
xlabel('Time(ms)');

stats = zeros(1,nwin);

for ipat = 1:np
    pat_list = 1:np;
    pat_list(ipat)=[];
    for iwin = 1:nwin
        [h,p, t, df] =  ttestch(mean(X2(pat_list,iwin)), std(X2(pat_list,iwin)), X2(ipat,iwin), np-1, 0.05);
        stats(iwin) = p;
        
    end
    stats = fdr(stats);
    fprintf('low_gt_high %d, numb of p<0.05 = %d\n',ipat, sum(stats<0.05));

end