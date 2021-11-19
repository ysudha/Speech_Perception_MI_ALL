% Same asR12, best adjacency matrices. Plotting them to visualize it
% (1) load graphs, edge = lowVOT ...looking at increases and decreases
% (2) load graphs, edge = midVOT ...looking at increases and decreases
% (3) load graphs, edge = highVOT ...looking at increases and decreases

 %% (1) load graphs, edge = lowVOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Robust PCA in time%%%%%%%%%%%%%%%%%%

miVOT123 = load('mi_matrix_VOT123_K4.mat'); %with mem 2
miVOT456 = load('mi_matrix_VOT456_K4.mat'); 
np = 22; %number of people 21 undergrads, 1 aphasia
nwin = 43;
nch = 30;
ntim = [1 51 101:5:301];
% creating the adjacency matrices, with edge = VOT456- VOT123, only
diff_mid_lowVOT = zeros(np,nch,nch,nwin-2);

for ipat = 1:np
    mi123 = miVOT123.mi_matrix_VOT123_K4{ipat,1}; %contains 43 time windows
    mi456 = miVOT456.mi_matrix_VOT456_K4{ipat,1};
    for iwin = 3:nwin
       mat1 = squeeze(mi456{iwin,1} - mi456{1,1}); 
       mat2 = squeeze(mi123{iwin,1} - mi123{1,1}); 
       diff_mid_lowVOT(ipat,:,:,iwin-2) = squeeze(mat1-mat2);
    end

end

nwin = nwin-2;
for ipat = 1:np
    figure; hold on;
    X = squeeze(diff_mid_lowVOT(ipat,:,:,:));
    for iwin = 1:nwin
        subplot(7,7,iwin);
        mat = squeeze(X(:,:,iwin));
        imagesc(mat);
        colormap jet; caxis([-0.3 0.3]);
        axis square;
    end
end



%% (2) load graphs, edge = midVOT - highVOT ...looking at increases and decreases
%%%% load new graphs, edge = midVOT - highVOT ...looking at increases only
miVOT789 = load('mi_matrix_VOT789_K4.mat'); %with mem 2
miVOT456 = load('mi_matrix_VOT456_K4.mat'); 
np = 22; %number of people 21 undergrads, 1 aphasia
nwin = 43;
nch = 30;
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
       diff_midVOT(ipat,:,:,iwin-2) = squeeze(mat1);
       diff_mid_highVOT(ipat,:,:,iwin-2) = squeeze(mat1-mat2);
    end
end

nwin = nwin-2;
energy_S_mid_high = zeros(np,nwin);
energy_S_mid = zeros(np,nwin);

parfor ipat = 1:np
    tic
    X = squeeze(diff_mid_highVOT(ipat,:,:,:));
    [Lx,Sx] = RobustPCA_time(X);
    Y = squeeze(diff_midVOT(ipat,:,:,:));
    [Ly,Sy] = RobustPCA_time(Y);
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
plot(X1(21,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S mid-high');

figure; 
 plot(X2(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X2(21,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S just mid');

[coeff, score, latent, tsquared, explained, mu] = pca(X1(:,1:end));%20
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('baseline norm Edge = mid-high');

[coeff, score, latent, tsquared, explained, mu] = pca(X2(:,1:end));%20
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('baseline norm Edge = mid');

%% (3) load graphs, edge = lowVOT - highVOT ...looking at increases and decreases
%%%% load new graphs, edge = lowVOT - highVOT ...looking at increases only
miVOT123 = load('mi_matrix_VOT123_K4.mat'); %with mem 2
miVOT789 = load('mi_matrix_VOT789_K4.mat'); 
np = 22; %number of people 21 undergrads, 1 aphasia
nwin = 43;
nch = 30;
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
       diff_highVOT(ipat,:,:,iwin-2)= squeeze(mat2);
       diff_low_highVOT(ipat,:,:,iwin-2) = squeeze(mat1-mat2);
    end

end

nwin = nwin-2;
energy_S_low_high = zeros(np,nwin);
energy_S_high = zeros(np,nwin);
S_pat = cell(1,np);
L_pat = cell(1,np);
parfor ipat = 1:np
    tic
    X = squeeze(diff_low_highVOT(ipat,:,:,:));
    [Lx,Sx] = RobustPCA_time(X);
    Y = squeeze(diff_highVOT(ipat,:,:,:));
    [Ly,Sy] = RobustPCA_time(Y);
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
plot(X1(21,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S low-high');

figure; 
 plot(X2(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X2(21,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S just high');

[coeff, score, latent, tsquared, explained, mu] = pca(X1(:,1:end));%20
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('baseline norm Edge = low-high');

[coeff, score, latent, tsquared, explained, mu] = pca(X2(:,1:end));%20
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('baseline norm Edge = high');

%% plotting
figure; hold on;
for ipat = 1:np

    X = squeeze(diff_mid_lowVOT(ipat,:,:,:));
    [Ly,Sy] = RobustPCA_time(X);
    [U2,S2,V2] = eig(Ly);
    sum(abs(diag(S2)))
    for iwin = 20%1:nwin
        subplot(5,5,ipat);
        mat = squeeze(Sy(:,:,iwin));
        [U1,S1,V1] = eig(mat);
        imagesc(U1(:,1:1)*S1(1:1,1:1)*V1(:,1:1)');
        colormap jet; caxis([0 0.7]);
        axis square;
        str = sprintf('%0.2f',energy_S_ml_time(ipat,iwin));
        title(str);
    end
end
