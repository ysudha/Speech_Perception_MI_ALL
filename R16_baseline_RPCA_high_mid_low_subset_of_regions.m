% do incerases in MI in (mid-high) or (mid-low) or (high-low) matter more or do
% decreases matter more?
%Just like R13, but looking at frontal+temporal Left+right
% (1) load graphs, edge = lowVOT ...looking at increases, decreases
% (2) load graphs, edge = midVOT ...looking at increases, decreases
% (3) load graphs, edge = highVOT ...looking at increases, decreases

 %% (1) load graphs, edge = lowVOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Robust PCA in time%%%%%%%%%%%%%%%%%%

miVOT123 = load('mi_matrix_VOT123_K4.mat'); %with mem 2
miVOT456 = load('mi_matrix_VOT456_K4.mat'); 
np = 22; %number of people 21 undergrads, 1 aphasia
nwin = 43;
%nch = 30;
ntim = [1 51 101:5:301];
% creating the adjacency matrices, with edge = VOT456- VOT123, only
% positive increases diff_mid_lowVOT

chlist = [1 3 4 8 12 13 2 6 7 11 15 16];
nch = length(chlist);

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
parfor ipat = 1:np
    tic
    X = squeeze(diff_mid_lowVOT(ipat,:,:,:));
    X(X<0)=0;
    [Lx,Sx] = RobustPCA_time(X);
    Y = squeeze(diff_mid_lowVOT(ipat,:,:,:));
    Y(Y>0)=0;
    [Ly,Sy] = RobustPCA_time(Y);
    for iwin = 1:nwin
        [U1,S1,V1] = eig(squeeze(Sx(:,:,iwin)));
        energy_S_mid_gt_low(ipat,iwin) = sum(abs(diag(S1)));
        [U2,S2,V2] = eig(squeeze(Sy(:,:,iwin)));
        energy_S_mid_less_low(ipat,iwin) = sum(abs(diag(S2)));
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

%% (2) load graphs, edge = midVOT - highVOT ...looking at increases and decreases
%%%% load new graphs, edge = midVOT - highVOT ...looking at increases only
miVOT789 = load('mi_matrix_VOT789_K4.mat'); %with mem 2
miVOT456 = load('mi_matrix_VOT456_K4.mat'); 
np = 22; %number of people 21 undergrads, 1 aphasia
nwin = 43;
nch = 30;
ntim = [1 51 101:5:301];
chlist = [1 3 4 8 12 13 2 6 7 11 15 16];
nch = length(chlist);
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

%% (3) load graphs, edge = lowVOT - highVOT ...looking at increases and decreases
%%%% load new graphs, edge = lowVOT - highVOT ...looking at increases only
miVOT123 = load('mi_matrix_VOT123_K4.mat'); %with mem 2
miVOT789 = load('mi_matrix_VOT789_K4.mat'); 
np = 22; %number of people end undergrads, 1 aphasia
nwin = 43;
nch = 30;
ntim = [1 51 101:5:301];
chlist = [1 3 4 8 12 13 2 6 7 11 15 16];
nch = length(chlist);
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
