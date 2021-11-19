% baseline norm and differences in VOT, decompose each matrix into L
% and S. Look at energy in both. - didnt seem as good as those in R5, that
% were NOT baseline normalized. 
% (1) load graphs, edge = midVOT - lowVOT ...looking at increases and decreases
% (2) load graphs, edge = midVOT - highVOT ...looking at increases and decreases
% (3) load graphs, edge = lowVOT - highVOT ...looking at increases and decreases

 %% (1) load graphs, edge = BSmidVOT - BSlowVOT ...looking at increases and decreases
%%%% load new graphs, edge = midVOT - lowVOT ...looking at increases only
miVOT123 = load('mi_matrix_VOT123_K4.mat'); %with mem 2
miVOT456 = load('mi_matrix_VOT456_K4.mat'); 
np = 22; %number of people 21 undergrads, 1 aphasia
nwin = 43;
nch = 30;
ntim = [1 51 101:5:301];
% creating the adjacency matrices, with edge = VOT456- VOT123, only
% positive increases diff_mid_lowVOT
diff_mid_lowVOT = zeros(np,nch,nch,nwin-2);
energy_L = zeros(np,nwin-2);
energy_S = zeros(np,nwin-2);
for ipat = 1:np
    mi123 = miVOT123.mi_matrix_VOT123_K4{ipat,1}; %contains 43 time windows
    mi456 = miVOT456.mi_matrix_VOT456_K4{ipat,1};
    tic
    parfor iwin = 3:nwin
        mat1 = squeeze(mi456{iwin,1} - mi456{1,1}); 
       mat2 = squeeze(mi123{iwin,1} - mi123{1,1}); 
       mat = squeeze(mat1-mat2);
       diff_mid_lowVOT(ipat,:,:,iwin-2)= mat;
       [L,S] = RobustPCA(mat);
       [U1,S1,V1] = eig(L);
       [U2,S2,V2] = eig(S);
       energy_L(ipat,iwin-2) = sum(abs(diag(S1)));
       energy_S(ipat,iwin-2) = sum(abs(diag(S2)));
    end
    toc
    
end

% plotting
X1 = energy_L(:,1:end);
X2 = energy_S(:,1:end);

figure; 
 plot(X1(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X1(end,:)','-r*','LineWidth',4);hold on;
title('Energy of L');

figure; 
 plot(X2(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X2(end,:)','-r*','LineWidth',4);hold on;
title('Energy of S');

[coeff, score, latent, tsquared, explained, mu] = pca(X1(:,1:end));%20
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
title('PCA of L');

[coeff, score, latent, tsquared, explained, mu] = pca(X2(:,1:end));%20
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
title('PCA of S');
% 
% 
% 
% %% (2) load graphs, edge = midVOT - highVOT ...looking at increases and decreases
% %%%% load new graphs, edge = midVOT - highVOT ...looking at increases only
% miVOT789 = load('mi_matrix_VOT789_K4.mat'); %with mem 2
% miVOT456 = load('mi_matrix_VOT456_K4.mat'); 
% np = 22; %number of people end undergrads, 1 aphasia
% nwin = 43;
% nch = 30;
% ntim = [1 51 101:5:301];
% % creating the adjacency matrices, with edge = VOT456- VOT123, only
% % positive increases diff_mid_lowVOT
% diff_mid_highVOT = zeros(np,nch,nch,nwin-2);
% energy_L = zeros(np,nwin);
% energy_S = zeros(np,nwin);
% for ipat = 1:np
%     mi789 = miVOT789.mi_matrix_VOT789_K4{ipat,1}; %contains 43 time windows
%     mi456 = miVOT456.mi_matrix_VOT456_K4{ipat,1};
%     tic
%     for iwin = 3:nwin
%         mat1 = squeeze(mi456{iwin,1} - mi456{1,1}); 
%        mat2 = squeeze(mi789{iwin,1} - mi789{1,1}); 
%        mat = squeeze(mat1-mat2);
%        diff_mid_highVOT(ipat,:,:,iwin-2)= mat;
% %        [L,S] = RobustPCA(mat);
% %        [U1,S1,V1] = eig(L);
% %        [U2,S2,V2] = eig(S);
% %        energy_L(ipat,iwin) = sum(abs(diag(S1)));
% %        energy_S(ipat,iwin) = sum(abs(diag(S2)));
%     end
%     toc
% end
% 
% % plotting
% X1 = energy_L(:,3:end);
% X2 = energy_S(:,3:end);
% 
% figure; 
%  plot(X1(1:end-1,:)','-b');hold on;
% plot(X1(end,:)','-r');hold on;
% title('Energy of L');
% 
% figure; 
%  plot(X2(1:end-1,:)','-b');hold on;
% plot(X2(end,:)','-r');hold on;
% title('Energy of S');
% 
% % %X1(10,:)=[];
% [coeff, score, latent, tsquared, explained, mu] = pca(X2(:,1:end));%20
% figure; scatter(score(1:end-1,1),score(1:end-1,2),'b*');hold on;scatter(score(end,1),score(end,2),'r*');
% figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%Robust PCA in time%%%%%%%%%%%%%%%%%%
% %%
% nwin = nwin-2;
% energy_Stime = zeros(np,nwin);
% card_Stime = zeros(np,nwin);
% parfor ipat = 1:np
%     tic
%     X = squeeze(diff_mid_highVOT(ipat,:,:,:));
%     [L,S] = RobustPCA_time(X);
%     for iwin = 1:nwin
%         [U2,S2,V2] = eig(squeeze(S(:,:,iwin)));
%         energy_Stime(ipat,iwin) = sum(abs(diag(S2)));
%         card_Stime(ipat,iwin) = nnz(squeeze(S(:,:,iwin)));
%     end
%     toc
% end
% 
% X1 = energy_Stime(:,3:end);
% X2 = card_Stime(:,3:end);
% E2 = X1;
% figure; 
%  plot(X1(1:end-1,:)','-b*','LineWidth',4);hold on;
% plot(X1(end,:)','-r*','LineWidth',4);hold on;
% title('Energy of S');
% 
% figure; 
%  plot(X2(1:end-1,:)','-b');hold on;
% plot(X2(end,:)','-r');hold on;
% title('Card of S');
% 
% [coeff, score, latent, tsquared, explained, mu] = pca(X1(:,1:end));%20
% figure; scatter(score(1:end-1,1),score(1:end-1,2),'b*');hold on;scatter(score(end,1),score(end,2),'r*');
% figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
% xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
% 
% 
% 
% 
% %% (3) load graphs, edge = lowVOT - highVOT ...looking at increases and decreases
% %%%% load new graphs, edge = lowVOT - highVOT ...looking at increases only
% miVOT123 = load('mi_matrix_VOT123_K4.mat'); %with mem 2
% miVOT789 = load('mi_matrix_VOT789_K4.mat'); 
% np = 22; %number of people end undergrads, 1 aphasia
% nwin = 43;
% nch = 30;
% ntim = [1 51 101:5:301];
% % creating the adjacency matrices, with edge = VOT456- VOT123, only
% % positive increases diff_mid_lowVOT
% diff_high_lowVOT = zeros(np,nch,nch,nwin-2);
% energy_L = zeros(np,nwin);
% energy_S = zeros(np,nwin);
% for ipat = 1:np
%     mi123 = miVOT123.mi_matrix_VOT123_K4{ipat,1}; %contains 43 time windows
%     mi789 = miVOT789.mi_matrix_VOT789_K4{ipat,1};
%     tic
%     for iwin = 3:nwin
%        mat1 = squeeze(mi123{iwin,1} - mi123{1,1}); 
%        mat2 = squeeze(mi789{iwin,1} - mi789{1,1}); 
%        mat = squeeze(mat2-mat1);
%        diff_high_lowVOT(ipat,:,:,iwin-2)= mat;
% 
% %        [L,S] = RobustPCA(mat);
% %        [U1,S1,V1] = eig(L);
% %        [U2,S2,V2] = eig(S);
% %        energy_L(ipat,iwin) = sum(abs(diag(S1)));
% %        energy_S(ipat,iwin) = sum(abs(diag(S2)));
%     end
%     toc
% end
% 
% 
% % plotting
% X1 = energy_L(:,3:end);
% X2 = energy_S(:,3:end);
% 
% figure; 
%  plot(X1(1:end-1,:)','-b');hold on;
% plot(X1(end,:)','-r');hold on;
% title('Energy of L');
% 
% figure; 
%  plot(X2(1:end-1,:)','-b');hold on;
% plot(X2(end,:)','-r');hold on;
% title('Energy of S');
% 
% % %X1(10,:)=[];
% [coeff, score, latent, tsquared, explained, mu] = pca(X2(:,1:end));%20
% figure; scatter(score(1:end-1,1),score(1:end-1,2),'b*');hold on;scatter(score(end,1),score(end,2),'r*');
% figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%Robust PCA in time%%%%%%%%%%%%%%%%%%
% %%
% nwin = nwin-2;
% energy_Stime = zeros(np,nwin);
% card_Stime = zeros(np,nwin);
% parfor ipat = 1:np
%     tic
%     X = squeeze(diff_high_lowVOT(ipat,:,:,:));
%     [L,S] = RobustPCA_time(X);
%     for iwin = 1:nwin
%         [U2,S2,V2] = eig(squeeze(S(:,:,iwin)));
%         energy_Stime(ipat,iwin) = sum(abs(diag(S2)));
%         card_Stime(ipat,iwin) = nnz(squeeze(S(:,:,iwin)));
%     end
%     toc
% end
% 
% X1 = energy_Stime(:,1:end);
% X2 = card_Stime(:,1:end);
% E3 = X1;
% figure; 
%  plot(X1(1:end-1,:)','-b*','LineWidth',4);hold on;
% plot(X1(end,:)','-r*','LineWidth',4);hold on;
% title('Energy of S');
% 
% figure; 
%  plot(X2(1:end-1,:)','-b');hold on;
% plot(X2(end,:)','-r');hold on;
% title('Card of S');
% 
% [coeff, score, latent, tsquared, explained, mu] = pca(X1(:,1:end));%20
% figure; scatter(score(1:end-1,1),score(1:end-1,2),'b*');hold on;scatter(score(end,1),score(end,2),'r*');
% figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
% xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
% 
% 
% 
% 
% 
% %% Combine the error vectors from all conditions
% 
% X = [E1 E2 E3];
% 
% [coeff, score, latent, tsquared, explained, mu] = pca(X(:,1:end));%20
% 
% figure; scatter(score(1:end-1,1),score(1:end-1,2),'b*');hold on;scatter(score(end,1),score(end,2),'r*');
% figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
% xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
% 
% Xhat = score(:,1:3)*coeff(:,1:3)' + mu;
%  figure; imagesc(X); figure; imagesc(Xhat)
%  
%  figure;
%  plot(Xhat(1:end-1,:)','-b*','LineWidth',4);hold on;
% plot(Xhat(end,:)','-r*','LineWidth',4);hold on;
% 
